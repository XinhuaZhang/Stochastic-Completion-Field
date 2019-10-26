-- This module has functions that store pinwheels to disk and read them out to do convolution,
-- for the purpose of handling the case when filter size is greater than RAM size.
module STC.Binary
  ( makePlanBinary
  , computeInitialEigenVectorBinary
  , writeDFTPinwheel
  , convolutionBinary
  ) where

import           Array.UnboxedArray           as AU
import           Control.Monad                as M
import           Control.Monad.IO.Class
import           Control.Monad.Trans.Resource
import           Data.Array.Repa              as R
import           Data.Binary
import           Data.ByteString.Lazy         as BSL
import           Data.ByteString.Lazy         as BS
import           Data.Complex
import           Data.Conduit
import           Control.DeepSeq              (($!!))
import           Data.Conduit.List            as CL
import           Data.List                    as L
import           Data.Vector.Storable         as VS
import           Data.Vector.Unboxed          as VU
import           DFT.Plan
import           Filter.Pinwheel
import           FokkerPlanck.Interpolation
import           FokkerPlanck.Pinwheel        (cutoff)
import           GHC.Float
import           System.IO
import           System.Random
import           Types
import           Utils.Array
import           Utils.Parallel               (ParallelParams (..),
                                               parConduitIO)

makePlanBinary ::
     DFTPlan -> Bool -> FilePath -> Int -> Int -> Int -> Int -> IO DFTPlan
makePlanBinary oldPlan wisdomFlag wisdomFilePath numThetaFreqs numScaleFreqs rows cols = do
  let n = numThetaFreqs * numScaleFreqs * rows * cols
  xs <- M.replicateM n randomIO :: IO [Double]
  ys <- M.replicateM n randomIO :: IO [Double]
  let vecTemp = VS.fromList . L.zipWith (:+) xs $ ys
  when wisdomFlag (importFFTWWisdom wisdomFilePath)
  lock <- getFFTWLock
  plan <-
    fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [cols, rows, numThetaFreqs, numScaleFreqs]
       [0, 1]
       vecTemp >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [cols, rows, numThetaFreqs, numScaleFreqs]
         [0, 1]
         vec >>= \(plan, vec) ->
         dft1dGPlan
           lock
           plan
           [numThetaFreqs, numScaleFreqs, cols, rows]
           [0, 1]
           vec >>= \(plan, vec) ->
           idft1dGPlan
             lock
             plan
             [numThetaFreqs, numScaleFreqs, cols, rows]
             [0, 1]
             vec)
  exportFFTWWisdom wisdomFilePath
  return plan


computeInitialEigenVectorBinary ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> [R2S1RPPoint]
  -> R2Z2Array
computeInitialEigenVectorBinary xLen yLen thetaFreqs scaleFreqs xs =
  let numThetaFreqs = L.length thetaFreqs
      numScaleFreqs = L.length scaleFreqs
      xShift = div xLen 2
      (xMin, xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      yShift = div yLen 2
      (yMin, yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      vec =
        toUnboxedVector .
        AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
        L.map
          (\(R2S1RPPoint (x, y, _, _)) ->
             ((x, y), 1 / (fromIntegral . L.length $ xs))) $
        xs
  in computeS .
     R.traverse
       (fromUnboxed (Z :. xLen :. yLen) vec)
       (const (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen)) $ \f idx@(Z :. tf :. sf :. i :. j) ->
       if tf == div numThetaFreqs 2 && sf == div numScaleFreqs 2
         then f (Z :. i :. j)
         else 0

{-# INLINE convertDFTPinwheel #-}
convertDFTPinwheel ::
     DFTPlan
  -> R.Array D DIM5 Double
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> Double
  -> (Int, Int)
  -> (Int, Int)
  -> IO BS.ByteString
convertDFTPinwheel plan radialArr thetaFreqs scaleFreqs hollowRadius rMax (rows, cols) (t, s) =
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. _) =
        extent radialArr
      pinwheelArr =
        traverse2
          thetaFreqs
          scaleFreqs
          (\_ _ -> (Z :. numTheta0Freq :. numScale0Freq :. cols :. rows)) $ \ft0 fs0 (Z :. t0 :. s0 :. i :. j) ->
          pinwheelFunc
            (PinwheelHollow0 hollowRadius)
            (thetaFreqs R.! (Z :. t) - ft0 (Z :. t0))
            (scaleFreqs R.! (Z :. s) + fs0 (Z :. s0))
            rMax
            0
            (i - center cols)
            (j - center rows)
      interpolatedPinwheelArr =
        radialCubicInterpolation
          (R.slice radialArr (Z :. t :. s :. All :. All :. All))
          1
          pinwheelArr
  in do xx <-
          fmap
            (encode .
             L.map (\(x :+ y) -> double2Float x :+ double2Float y) . VS.toList) .
          dftExecute
            plan
            (DFTPlanID DFT1DG [cols, rows, numTheta0Freq, numScale0Freq] [0, 1]) .
          VS.convert . toUnboxed . computeS . rotate4D2 . makeFilter2D $
          interpolatedPinwheelArr
        return $!! xx

{-# INLINE writeDFTPinwheelSink #-}
writeDFTPinwheelSink ::
     FilePath -> ConduitT BS.ByteString Void (ResourceT IO) ()
writeDFTPinwheelSink filePath = do
  h <- liftIO $ openBinaryFile filePath WriteMode
  go h
  where
    go h = do
      z <- await
      case z of
        Nothing -> liftIO $ hClose h
        Just bs -> do
          liftIO . BSL.hPut h . encode . BS.length $ bs
          liftIO . BS.hPut h $ bs
          go h

writeDFTPinwheel ::
     ParallelParams
  -> DFTPlan
  -> R.Array D DIM5 Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> Double
  -> (Int, Int)
  -> FilePath
  -> IO ()
writeDFTPinwheel parallelParams plan radialArr hollowRadius cutoffRadius thetaFreqs scaleFreqs rMax (rows, cols) filePath = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. _) =
        extent radialArr
      thetaFreqsArr = fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs
      scaleFreqsArr = fromListUnboxed (Z :. L.length scaleFreqs) scaleFreqs
  runConduitRes $
    CL.sourceList
      [(t, s) | t <- [0 .. numThetaFreq - 1], s <- [0 .. numScaleFreq - 1]] .|
    parConduitIO
      parallelParams
      (convertDFTPinwheel
         plan
         (cutoff (round cutoffRadius) radialArr)
         thetaFreqsArr
         scaleFreqsArr
         hollowRadius
         rMax
         (rows, cols)) .|
    writeDFTPinwheelSink filePath

{-# INLINE sourceBinary #-}
sourceBinary :: FilePath -> ConduitT () BS.ByteString (ResourceT IO) ()
sourceBinary filePath = do
  h <- liftIO $ openBinaryFile filePath ReadMode
  go h
  where
    go handle = do
      lenBS <- liftIO (BSL.hGet handle 8)
      if BSL.null lenBS
        then liftIO $ hClose handle
        else do
          let len = decode lenBS :: Int
          bs <- liftIO . BS.hGet handle $ len
          yield bs
          go handle

{-# INLINE convolution #-}
convolution ::
     Bool
  -> DFTPlan
  -> DIM4
  -> R.Array U DIM1 Double
  -> VS.Vector (Complex Double)
  -> BS.ByteString
  -> IO (VU.Vector (Complex Double))
convolution True plan dim@(Z :. cols :. rows :. numThetaFreq :. numScaleFreq) thetaFreqs input bs = do
  arr <-
    fmap (fromUnboxed dim . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [cols, rows, numThetaFreq, numScaleFreq] [0, 1]) .
    VS.zipWith (*) input $
    VS.fromList . L.map (\(x :+ y) -> float2Double x :+ float2Double y) . decode $
    bs
  return $!! toUnboxed . sumS . sumS . R.traverse2 arr thetaFreqs const $
    (\f ft idx@(Z :. _ :. _ :. i :. _) -> f idx * exp (0 :+ ft (Z :. i) * pi))
convolution False plan dim@(Z :. cols :. rows :. numThetaFreq :. numScaleFreq) thetaFreqs input bs = do
  x <-
    fmap (toUnboxed . sumS . sumS . fromUnboxed dim . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [cols, rows, numThetaFreq, numScaleFreq] [0, 1]) .
    VS.zipWith (*) input .
    VS.fromList . L.map (\(x :+ y) -> float2Double x :+ float2Double y) . decode $
    bs
  return $!! x


convolutionBinary ::
     Bool
  -> ParallelParams
  -> DFTPlan
  -> FilePath
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
convolutionBinary sinkFlag parallelParams plan filterFilePath thetaFreqs inputArr = do
  let thetaFreqsArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
      (Z :. numThetaFreq :. numScaleFreq :. cols :. rows) = extent inputArr
  inputArrF <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [cols, rows, numThetaFreq, numScaleFreq] [0, 1]) .
    VU.convert . toUnboxed . computeS . rotate4D2 $
    inputArr
  xs <-
    runConduitRes $
    sourceBinary filterFilePath .|
    parConduitIO
      parallelParams
      (convolution
         sinkFlag
         plan
         (Z :. cols :. rows :. numThetaFreq :. numScaleFreq)
         thetaFreqsArr
         inputArrF) .|
    CL.consume
  let arr = fromUnboxed (extent inputArr) . VU.concat $ xs
  return $
    if sinkFlag
      then computeS . R.traverse2 arr thetaFreqsArr const $ \f ft idx@(Z :. i :. _ :. _ :. _) ->
             f idx * exp (0 :+ ft (Z :. i) * pi)
      else arr
