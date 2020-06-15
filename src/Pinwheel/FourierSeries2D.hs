{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Strict              #-}
module Pinwheel.FourierSeries2D where

import           Control.Monad                  as M
import           Control.Monad.IO.Class
import           Control.Monad.Trans.Resource
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX as A
import           Data.Array.IArray              as IA
import           Data.Array.Repa                as R
import           Data.Complex
import           Data.Conduit                   as C
import           Data.Conduit.List              as CL
import           Data.List                      as L
import           Data.Vector.Generic            as VG
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.FourierSeries2D
import           Math.Gamma
import           Pinwheel.Base
import           Pinwheel.List
import           Utils.BLAS
import           Utils.Distribution
import           Utils.List
import           Utils.Parallel                 hiding ((.|))
import           Utils.SimpsonRule
import           Utils.Time

pinwheelFourierSeries ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> Int
  -> Int
  -> IO (IA.Array (Int, Int) (vector (Complex e)))
pinwheelFourierSeries deviceIDs ptxs numR2Freqs numPoints delta periodR2 phiFreq rhoFreq thetaFreq rFreq sigma numBatchR2 numBatchPinwheelFreqs = do
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- pinwheelFreqs rhoFreq rFreq
        , angularFreq <- pinwheelFreqs phiFreq thetaFreq
        ]
      centerR2Freq = div numR2Freqs 2
  pinwheels <-
    computeUnboxedP .
    R.traverse
      (fromListUnboxed (Z :. (L.length idxs)) idxs)
      (\(Z :. freqs) -> (Z :. freqs :. numR2Freqs :. numR2Freqs)) $ \f (Z :. pFreq :. xFreq :. yFreq) ->
      fourierMellin
        sigma
        (snd $ f (Z :. pFreq))
        (fst $ f (Z :. pFreq))
        ( delta * fromIntegral (xFreq - centerR2Freq)
        , delta * fromIntegral (yFreq - centerR2Freq))
  -- pinwheels <-
  --   M.mapM
  --     (\idx ->
  --        fmap
  --          (CuMat (numR2Freqs ^ 2) (L.length idx) .
  --           CuVecHost . VU.convert . toUnboxed) .
  --        computeP .
  --        R.traverse
  --          (fromListUnboxed (Z :. (L.length idx)) idx)
  --          (\(Z :. freqs) -> (Z :. numR2Freqs :. numR2Freqs :. freqs)) $ \f (Z :. xFreq :. yFreq :. pFreq) ->
  --          fourierMellin
  --            sigma
  --            (snd $ f (Z :. pFreq))
  --            (fst $ f (Z :. pFreq))
  --            ( fromIntegral (xFreq - centerR2Freq)
  --            , fromIntegral (yFreq - centerR2Freq))) .
  --   divideListN numBatchPinwheelFreqs $
  --   idxs
  -- printCurrentTime "Compute Fourier Series... "
  -- arr <-
  --   computeFourierSeriesR2Stream
  --     deviceIDs
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     delta
  --     numBatchR2
  --     pinwheels
  -- printCurrentTime "Compute Fourier Series Done"
  let (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound phiFreq thetaFreq
  return .
    listArray ((radialLB, angularLB), (radialUB, angularUB)) .
    parMap
      rdeepseq
      (\i ->
         VG.convert . toUnboxed . computeS . R.slice pinwheels $  -- arr $
         (Z :. i :. All :. All)) $
    [0 .. (L.length idxs - 1)]


-- Stream
{-# INLINE source #-}
source :: Int -> Int -> Int -> Int -> ConduitT () (Int, Int) (ResourceT IO) ()
source phiFreq rhoFreq thetaFreq rFreq =
  CL.sourceList
    [ (radialFreq, angularFreq)
    | radialFreq <- pinwheelFreqs rhoFreq rFreq
    , angularFreq <- pinwheelFreqs phiFreq thetaFreq
    ]

conduit ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , VG.Vector vector (Complex e)
     )
  => [Int]
  -> Int
  -> Int
  -> e
  -> Int
  -> e
  -> [[CuMat (Complex e)]]
  -> ConduitT (Int, Int) [vector (Complex e)] (ResourceT IO) ()
conduit deviceIDs numR2Freqs numPoints periodR2 batchSize sigma inverseR2Harmonics = do
  liftIO $ printCurrentTime "conduit"
  idx <- CL.take batchSize
  unless
    (L.null idx)
    (do let centerR2Freq = div numR2Freqs 2
        pinwheels <-
          fmap
            (CuMat (L.length idx) (numR2Freqs ^ 2) .
             CuVecHost . VU.convert . toUnboxed) .
          liftIO .
          computeP .
          R.traverse
            (fromListUnboxed (Z :. (L.length idx)) idx)
            (\(Z :. freqs) -> (Z :. freqs :. numR2Freqs :. numR2Freqs)) $ \f (Z :. pFreq :. xFreq :. yFreq) ->
            fourierMellin
              sigma
              (snd $ f (Z :. pFreq))
              (fst $ f (Z :. pFreq))
              ( fromIntegral (xFreq - centerR2Freq)
              , fromIntegral (yFreq - centerR2Freq))
        arr <-
          liftIO $
          computeFourierSeriesR2
            deviceIDs
            numR2Freqs
            numPoints
            periodR2
            inverseR2Harmonics
            [pinwheels]
        yield .
          L.map
            (\i ->
               VG.convert . toUnboxed . computeS . R.slice arr $
               (Z :. i :. All :. All)) $
          [0 .. (L.length idx - 1)])

pinwheelFourierSeriesStream ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , VG.Vector vector (Complex e)
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> Int
  -> Int
  -> Int
  -> IO (IA.Array (Int, Int) (vector (Complex e)))
pinwheelFourierSeriesStream deviceIDs ptxs numR2Freqs numPoints delta periodR2 phiFreq rhoFreq thetaFreq rFreq sigma numBatchR2 batchSize numBatchPinwheelFreqs = do
  inverseR2Harmonics <-
    createInverseHarmonicMatriesGPU
      ptxs
      numBatchR2
      numPoints
      numR2Freqs
      periodR2
      delta
  let (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound rhoFreq rFreq
  fmap (listArray ((radialLB, angularLB), (radialUB, radialUB)) . L.concat) .
    runConduitRes $
    source phiFreq rhoFreq thetaFreq rFreq .|
    conduit
      deviceIDs
      numR2Freqs
      numPoints
      periodR2
      batchSize
      sigma
      inverseR2Harmonics .|
    CL.consume

-- 2D spatial domain to frequency domain:
-- transform: cis (-freq * x)   
-- inverse transform: cis (freq * x)   

-- Analytical solution type 1
-- The pinwheel harmonics are cis (-freq * x)   
{-# INLINE analyticalFourierSeriesFunc1 #-}
analyticalFourierSeriesFunc1 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e))
  => Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> Complex e
analyticalFourierSeriesFunc1 angularFreq radialFreq sigma period phi rho =
  ((0 :+ 1) ^ (abs angularFreq)) * ((cis (fromIntegral (-angularFreq) * phi))) *
  (period / rho :+ 0) *
  ((pi * rho / period :+ 0) ** ((-1 - sigma) :+ fromIntegral radialFreq)) *
  (gamma $
   ((2 + fromIntegral (abs angularFreq) + sigma) :+ fromIntegral (-radialFreq)) / 2) /
  (gamma $ ((fromIntegral (abs angularFreq) - sigma) :+ (fromIntegral radialFreq)) / 2) 

-- The pinwheel is in the frequency domain, the function computes the Fourier series
{-# INLINE analyticalFourierSeries1 #-}
analyticalFourierSeries1 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierSeries1 numPoints delta angularFreq radialFreq sigma period =
  let center = div numPoints 2
  in fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
       let x = fromIntegral $ i - center
           y = fromIntegral $ j - center
           rho = sqrt $ x ^ 2 + y ^ 2
           phi = atan2 y x
       in if rho <= 0
            then 0
            else analyticalFourierSeriesFunc1
                   angularFreq
                   radialFreq
                   sigma
                   period
                   phi
                   (rho * delta)

-- The pinwheel is in the spatial domain, the function computes the Fourier coefficients
{-# INLINE analyticalFourierCoefficients1 #-}
analyticalFourierCoefficients1 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierCoefficients1 numFreqs delta angularFreq radialFreq sigma period =
  let c = (-1) ^ (abs angularFreq) :+ 0
  in R.map (* c) $
     analyticalFourierSeries1 numFreqs delta angularFreq radialFreq sigma period


-- Analytical solution type 2
-- The pinwheel harmonics are cis (freq * x)   
{-# INLINE analyticalFourierSeriesFunc2 #-}
analyticalFourierSeriesFunc2 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e))
  => Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> Complex e
analyticalFourierSeriesFunc2 angularFreq radialFreq sigma period phi rho =
  ((0 :+ 1) ^ (abs angularFreq)) * ((cis (fromIntegral angularFreq * phi))) *
  (period / rho :+ 0) *
  ((pi * rho / period :+ 0) ** ((-1 - sigma) :+ fromIntegral (-radialFreq))) *
  (gamma $
   ((2 + fromIntegral (abs angularFreq) + sigma) :+ fromIntegral radialFreq) / 2) /
  (gamma $
   ((fromIntegral (abs angularFreq) - sigma) :+ (fromIntegral (-radialFreq))) /
   2) 

-- The pinwheel is in the frequency domain, the function computes the Fourier series
{-# INLINE analyticalFourierSeries2 #-}
analyticalFourierSeries2 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierSeries2 numPoints delta angularFreq radialFreq sigma period =
  let center = div numPoints 2
  in fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
       let x = fromIntegral $ i - center
           y = fromIntegral $ j - center
           rho = sqrt $ x ^ 2 + y ^ 2
           phi = atan2 y x
       in if rho <= 0
            then 0
            else analyticalFourierSeriesFunc2
                   angularFreq
                   radialFreq
                   sigma
                   period
                   phi
                   (rho * delta)

-- The pinwheel is in the spatial domain, the function computes the Fourier coefficients
{-# INLINE analyticalFourierCoefficients2 #-}
analyticalFourierCoefficients2 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierCoefficients2 numFreqs delta angularFreq radialFreq sigma period =
  let c = (-1) ^ (abs angularFreq) :+ 0
  in R.map (* c) $
     analyticalFourierSeries2 numFreqs delta angularFreq radialFreq sigma period


pinwheelFourierCoefficientsAnatical ::
     ( Unbox e
     , RealFloat e
     , Gamma (Complex e)
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> IA.Array (Int, Int) (vector (Complex e))
pinwheelFourierCoefficientsAnatical numR2Freqs phiFreq rhoFreq thetaFreq rFreq sigma periodR2 =
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- pinwheelFreqs rhoFreq rFreq
        , angularFreq <- pinwheelFreqs phiFreq thetaFreq
        ]
      centerR2Freq = div numR2Freqs 2
      (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound phiFreq thetaFreq
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             VG.convert . toUnboxed . computeS $
             analyticalFourierCoefficients2
               numR2Freqs
               1
               (angularFreq)
               (radialFreq)
               sigma
               periodR2)
          idxs
  in listArray ((radialLB, angularLB), (radialUB, angularUB)) pinwheels
