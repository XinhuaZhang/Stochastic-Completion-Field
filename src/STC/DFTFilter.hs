{-# LANGUAGE FlexibleContexts #-}
module STC.DFTFilter where

import           Control.Monad                as M
import           Control.Monad.IO.Class
import           Control.Monad.Parallel       as MP
import           Control.Monad.Trans.Resource
import           Data.Array.Repa              as R
import           Data.Binary
import           Data.ByteString.Lazy         as BL
import           Data.Complex
import           Data.Conduit
import           Data.Conduit.List            as C
import           Data.Int
import           Data.List                    as L
import           Data.Vector.Storable         as VS
import           Data.Vector.Unboxed          as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.Interpolation
import           System.IO
import           Utils.Array
import           Utils.Parallel               hiding ((.|))

{-# INLINE sourceFile #-}
sourceFile :: FilePath -> ConduitT () (VS.Vector (Complex Double)) (ResourceT IO) ()
sourceFile filePath = do
  h <- liftIO $ openBinaryFile filePath ReadMode
  sourceFunc h
  liftIO $ hClose h
  where
    sourceFunc handle = do
      flag <- liftIO $ hIsEOF handle
      if flag
        then sourceNull
        else do
          lenBS <- liftIO $ BL.hGet handle 8
          let len = fromIntegral (decode lenBS :: Int64) :: Int
          dataBS <- liftIO $ BL.hGet handle len
          yield . VS.fromList $ (decode dataBS :: [Complex Double])
          sourceFunc handle 

{-# INLINE sinkFile #-}
sinkFile :: FilePath -> ConduitT ByteString Void (ResourceT IO) ()
sinkFile filePath = do
  h <- liftIO $ openBinaryFile filePath WriteMode
  sinkFunc h
  where
    sinkFunc handle = do
      x <- await
      case x of
        Nothing -> liftIO $ hClose handle
        Just y -> do
          let len = BL.length y
          liftIO $ BL.hPut handle (encode len)
          liftIO $ BL.hPut handle y
          sinkFunc handle            

{-# INLINE interpolation #-}
interpolation ::
     (R.Source r Double)
  => DFTPlan
  -> (Double -> Double -> Double -> Double -> Int -> Int -> Complex Double)
  -> R.Array r DIM5 Double
  -> Int
  -> Int
  -> Double
  -> Double
  -> (VU.Vector Double)
  -> (VU.Vector Double)
  -> (VU.Vector Double)
  -> (VU.Vector Double)
  -> (Int, Int)
  -> IO ByteString
interpolation plan pinwheelFunc radialArr xLen yLen scaleFactor rMax thetaFreqs scaleFreqs theta0Freqs scale0Freqs (t, s) = do
  let tf = thetaFreqs VU.! t
      sf = scaleFreqs VU.! s
      pinwheelArr =
        traverse2
          (fromUnboxed (Z :. VU.length theta0Freqs) theta0Freqs)
          (fromUnboxed (Z :. VU.length scale0Freqs) scale0Freqs)
          (\(Z :. numTheta0Freq) (Z :. numScale0Freq) ->
             (Z :. numTheta0Freq :. numScale0Freq :. xLen :. yLen)) $ \ft0 fs0 (Z :. t0 :. s0 :. i :. j) ->
          pinwheelFunc
            (tf - ft0 (Z :. t0))
            (sf + fs0 (Z :. s0))
            rMax
            0
            (i - center xLen)
            (j - center yLen)
      filter =
        radialCubicInterpolation
          (R.slice radialArr (Z :. t :. s :. All :. All :. All))
          scaleFactor
          pinwheelArr
  filterF <-
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [VU.length theta0Freqs, VU.length scale0Freqs, xLen, yLen]
         [2, 3]) .
    VU.convert . toUnboxed . computeS $
    filter
  return . encode . VS.toList $ filterF

{-# INLINE dftFilter2File #-}
dftFilter2File ::
     (R.Source r Double)
  => ParallelParams
  -> DFTPlan
  -> (Double -> Double -> Double -> Double -> Int -> Int -> Complex Double)
  -> R.Array r DIM5 Double
  -> Int
  -> Int
  -> Double
  -> Double
  -> (VU.Vector Double)
  -> (VU.Vector Double)
  -> (VU.Vector Double)
  -> (VU.Vector Double)
  -> FilePath
  -> IO ()
dftFilter2File parallelParams plan pinwheelFunc radialArr xLen yLen scaleFactor rMax thetaFreqs scaleFreqs theta0Freqs scale0Freqs filePath =
  runConduitRes $
  sourceList
    [ (t, s)
    | t <- [0 .. VU.length thetaFreqs]
    , s <- [0 .. VU.length scaleFreqs]
    ] .|
  parConduitIO
    parallelParams
    (interpolation
       plan
       pinwheelFunc
       radialArr
       xLen
       yLen
       scaleFactor
       rMax
       thetaFreqs
       scaleFreqs
       theta0Freqs
       scale0Freqs) .|
  sinkFile filePath

{-# INLINE convolutionHelper #-}
convolutionHelper ::
     DFTPlan
  -> Int
  -> Int
  -> Int
  -> Int
  -> VS.Vector (Complex Double)
  -> VS.Vector (Complex Double)
  -> IO (VU.Vector (Complex Double))
convolutionHelper plan numThetaFreq numScaleFreq xLen yLen vecF1 vecF2 =
  toUnboxed .
  sumS .
  sumS .
  rotate4D .
  rotate4D .
  fromUnboxed (Z :. numThetaFreq :. numScaleFreq :. xLen :. yLen) . VS.convert <$>
  dftExecute
    plan
    (DFTPlanID IDFT1DG [numThetaFreq, numScaleFreq, xLen, yLen] [2, 3])
    (VS.zipWith (*) vecF1 vecF2)
   

{-# INLINE convolveSlowly #-}
convolveSlowly ::
     ParallelParams
  -> FilePath
  -> DFTPlan
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
convolveSlowly parallelParams filePath plan arr = do
  let (Z :. numThetaFreq :. numScaleFreq :. xLen :. yLen) = extent arr
  xs <-
    runConduitRes $
    sourceFile filePath .|
    parConduitIO
      parallelParams
      (convolutionHelper
         plan
         numThetaFreq
         numScaleFreq
         xLen
         yLen
         (VU.convert . toUnboxed $ arr)) .|
    C.consume
  return . fromUnboxed (extent arr) . VU.concat $ xs


{-# INLINE convolutionHelperSink #-}
convolutionHelperSink ::
     DFTPlan
  -> Int
  -> Int
  -> Int
  -> Int
  -> VU.Vector Double
  -> VS.Vector (Complex Double)
  -> ((Int, Int), VS.Vector (Complex Double))
  -> IO (VU.Vector (Complex Double))
convolutionHelperSink plan numThetaFreq numScaleFreq xLen yLen thetaFreqs vecF1 ((t, _), vecF2) = do
  let tf = thetaFreqs VU.! t
  arr <-
    fromUnboxed (Z :. numThetaFreq :. numScaleFreq :. xLen :. yLen) . VS.convert <$>
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreq, numScaleFreq, xLen, yLen] [2, 3])
      (VS.zipWith (*) vecF1 vecF2)
  return . toUnboxed . sumS . sumS . rotate4D . rotate4D . R.traverse arr id $ \f idx@(Z :. t0 :. _ :. i :. j) ->
    f idx * exp (0 :+ ((thetaFreqs VU.! t0) + tf) * pi)


{-# INLINE convolveSlowlySink #-}
convolveSlowlySink ::
     ParallelParams
  -> FilePath
  -> DFTPlan
  -> VU.Vector Double
  -> VU.Vector Double
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
convolveSlowlySink parallelParams filePath plan thetaFreqs scaleFreqs arr = do
  let (Z :. numThetaFreq :. numScaleFreq :. xLen :. yLen) = extent arr
  xs <-
    runConduitRes $
    sourceFile filePath .|
    mergeSource
      (sourceList
         [ (t, s)
         | t <- [0 .. VU.length thetaFreqs]
         , s <- [0 .. VU.length scaleFreqs]
         ]) .|
    parConduitIO
      parallelParams
      (convolutionHelperSink
         plan
         numThetaFreq
         numScaleFreq
         xLen
         yLen
         thetaFreqs
         (VU.convert . toUnboxed $ arr)) .|
    C.consume
  return . fromUnboxed (extent arr) . VU.concat $ xs
