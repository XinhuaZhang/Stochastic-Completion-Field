{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
module FourierPinwheel.AsteriskGaussian where

import           Data.Array.Repa             as R
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Generic         as VG
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           DFT.Plan
import           Filter.Utils
import           FourierPinwheel.Array
import           FourierPinwheel.Hypergeo1F1
import           Math.Gamma
import           Pinwheel.FourierSeries2D
import           Utils.Distribution
import           Utils.Parallel


-- The envelope in R2 is r ^ alpha, where -2 < alpha < -0.5
{-# INLINE asteriskGaussian #-}
asteriskGaussian ::
     ( VG.Vector vector (Complex e)
     , RealFloat e
     , NFData (vector (Complex e))
     , Gamma (Complex e)
     , Unbox e
     )
  => Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> [vector (Complex e)]
asteriskGaussian numR2Freqs thetaFreq alpha periodR2 periodEnv std =
  parMap
    rdeepseq
    (\angularFreq ->
       let arr =
             R.map
               (* (sqrt (pi / std) *
                   exp (fromIntegral angularFreq ^ 2 / (-4) / std) :+
                   0)) .
             centerHollowArray numR2Freqs $
             analyticalFourierCoefficients1
               numR2Freqs
               1
               angularFreq
               0
               alpha
               periodR2
               periodEnv
        in VG.convert . toUnboxed . computeS $ arr)
    [-thetaFreq .. thetaFreq]

{-# INLINE asteriskGaussianEnvelope #-}
asteriskGaussianEnvelope ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> R.Array s DIM2 (Complex Double)
  -> IO (VS.Vector (Complex Double))
asteriskGaussianEnvelope plan numR2Freqs thetaFreq alpha periodR2 periodEnv stdTheta filter = do
  let planID = DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]
      inversePlanID = DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1]
      vecs =
        asteriskGaussian numR2Freqs thetaFreq alpha periodR2 periodEnv stdTheta
  filterF <-
    dftExecute plan planID . VU.convert . toUnboxed . computeS . makeFilter2D $
    filter
  asteriskGaussianF <- dftExecuteBatchP plan planID vecs
  fmap VS.concat .
    dftExecuteBatchP plan inversePlanID .
    parMap rdeepseq (VS.zipWith (*) filterF) $
    asteriskGaussianF

-- The envelope is r^alpha X Gaussian
asteriskGaussian2 ::
     DFTPlan
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> IO (VS.Vector (Complex Double))
asteriskGaussian2 plan numR2Freqs thetaFreq alpha periodR2 periodEnv stdR2 stdTheta = do
  let centerFreq = div numR2Freqs 2
      gaussian2D =
        fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
          let i = i' - centerFreq
              j = j' - centerFreq
           in exp
                (pi * fromIntegral (i ^ 2 + j ^ 2) /
                 ((-1) * periodR2 ^ 2 * stdR2 ^ 2)) /
              (2 * pi * stdR2 ^ 2) :+
              0
  asteriskGaussianEnvelope
    plan
    numR2Freqs
    thetaFreq
    alpha
    periodR2
    periodEnv
    stdTheta
    gaussian2D

-- The envelope is r^alpha X LowPass
asteriskGaussianLowPass ::
     DFTPlan
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> IO (VS.Vector (Complex Double))
asteriskGaussianLowPass plan numR2Freqs thetaFreq alpha periodR2 periodEnv stdTheta radius = do
  let lpf =
        fromUnboxed (Z :. numR2Freqs :. numR2Freqs) $
        idealLowPassFilter radius periodR2 numR2Freqs
  asteriskGaussianEnvelope
    plan
    numR2Freqs
    thetaFreq
    alpha
    periodR2
    periodEnv
    stdTheta
    lpf


{-# INLINE asteriskGaussianFull #-}
asteriskGaussianFull ::
     ( VG.Vector vector (Complex e)
     , RealFloat e
     , NFData (vector (Complex e))
     , Gamma (Complex e)
     , Unbox e
     )
  => Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> e
  -> [vector (Complex e)]
asteriskGaussianFull numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv stdTheta stdR =
  let zeroVec = VG.replicate (numR2Freqs ^ 2) 0
   in parMap
        rdeepseq
        (\(radialFreq, angularFreq) ->
           let pinwheel =
                 analyticalFourierCoefficients1
                   numR2Freqs
                   1
                   angularFreq
                   radialFreq
                   alpha
                   periodR2
                   periodEnv
               -- arr =
               --   R.map
               --     (* ((-- gaussian1DFreq (fromIntegral angularFreq) stdTheta *
               --          gaussian1DFourierCoefficients
               --            (fromIntegral radialFreq)
               --            (log periodEnv)
               --            stdR
               --          :+
               --          0) -- *
               --         -- cis
               --         --   (2 * pi / log periodEnv * fromIntegral radialFreq *
               --         --    log 0.5)
               --        )) $
               --   -- if radialFreq == 0 && angularFreq == 0
               --   --   then pinwheel
               --   --   else
                   -- centerHollowArray numR2Freqs pinwheel
               arr =
                 R.map
                   (* ((1 / (1 + (2 * pi / log periodEnv * fromIntegral radialFreq )^2)) :+ 0)
                      ) $
                   centerHollowArray numR2Freqs pinwheel
            in if angularFreq == 0
                 then VG.convert . toUnboxed . computeS $ arr
                 else zeroVec)
        [ (radialFreq, angularFreq)
        | radialFreq <-  [-rFreq .. rFreq]
        , angularFreq <- [-thetaFreq .. thetaFreq]
        ]


{-# INLINE gaussianFull #-}
gaussianFull ::
     ( VG.Vector vector (Complex e)
     , RealFloat e
     , NFData (vector (Complex e))
     , Gamma (Complex e)
     , Unbox e
     )
  => Int
  -> Int
  -> Int
  -> e
  -> e
  -> [vector (Complex e)]
gaussianFull numR2Freqs thetaFreq rFreq periodR2 stdR2 =
  let zeroVec = VG.replicate (numR2Freqs ^ 2) 0
      periodEnv = periodR2 ^ 2 / 4
   in parMap
        rdeepseq
        (\(radialFreq, angularFreq) ->
           if angularFreq == 0 -- && radialFreq == 0
             then VG.convert . toUnboxed . computeS $
                  fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
                    let i = fromIntegral $ i' - div numR2Freqs 2
                        j = fromIntegral $ j' - div numR2Freqs 2
                     in (gaussian2DFourierCoefficients i j periodR2 stdR2 :+ 0) /
                        (1 :+ (2 * pi / log periodEnv * fromIntegral radialFreq)^2)
             else zeroVec)
        [ (radialFreq, angularFreq)
        | radialFreq <- [-rFreq .. rFreq]
        , angularFreq <- [-thetaFreq .. thetaFreq]
        ]


{-# INLINE asteriskGaussianFullEnvelope #-}
asteriskGaussianFullEnvelope ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> R.Array s DIM2 (Complex Double)
  -> IO (VS.Vector (Complex Double))
asteriskGaussianFullEnvelope plan numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv stdTheta stdR filter = do
  let planID = DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]
      inversePlanID = DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1]
      vecs =
        asteriskGaussianFull
          numR2Freqs
          thetaFreq
          rFreq
          alpha
          periodR2
          periodEnv
          stdTheta
          stdR
      -- vecs =
      --  gaussianFull
      --      numR2Freqs
      --      thetaFreq
      --      rFreq
      --      periodR2
      --      stdR
  filterF <-
    dftExecute plan planID . VU.convert . toUnboxed . computeS . makeFilter2D $
    filter
  asteriskGaussianF <- dftExecuteBatchP plan planID vecs
  fmap VS.concat .
    dftExecuteBatchP plan inversePlanID .
    parMap rdeepseq (VS.zipWith (*) filterF) $
    asteriskGaussianF
  -- return . VS.concat $ vecs


-- The envelope is r^alpha X LowPass
asteriskGaussianFullLowPass ::
     DFTPlan
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> IO (VS.Vector (Complex Double))
asteriskGaussianFullLowPass plan numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv stdTheta stdR radius = do
  let lpf =
        fromUnboxed (Z :. numR2Freqs :. numR2Freqs) $
        idealLowPassFilter radius periodR2 numR2Freqs
  asteriskGaussianFullEnvelope
    plan
    numR2Freqs
    thetaFreq
    rFreq
    alpha
    periodR2
    periodEnv
    stdTheta
    stdR
    lpf


-- The envelope is r^alpha X LowPass
asteriskGaussian2Full ::
     DFTPlan
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> IO (VS.Vector (Complex Double))
asteriskGaussian2Full plan numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv stdTheta stdR stdR2 = do
  let centerFreq = div numR2Freqs 2
      gaussian2D =
        fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
          let i = fromIntegral $ i' - centerFreq
              j = fromIntegral $ j' - centerFreq
           in gaussian2DFourierCoefficients i j periodR2 stdR2 :+ 0
  asteriskGaussianFullEnvelope
    plan
    numR2Freqs
    thetaFreq
    rFreq
    alpha
    periodR2
    periodEnv
    stdTheta
    stdR
    gaussian2D


asteriskGaussianRFull ::
     DFTPlan
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> IO (VS.Vector (Complex Double))
asteriskGaussianRFull plan numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv stdTheta stdR stdR2 = do
  let planID = DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]
      inversePlanID = DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1]
      centerFreq = div numR2Freqs 2
      gaussian2D =
        fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
          let i = i' - centerFreq
              j = j' - centerFreq
           in exp
                (pi * fromIntegral (i ^ 2 + j ^ 2) /
                 ((-1) * periodR2 ^ 2 * stdR2 ^ 2)) /
              (2 * pi * stdR2 ^ 2) :+
              0
      r =
        centerHollowArray' numR2Freqs $ analyticalFourierCoefficients3 numR2Freqs 1 0 0 alpha periodR2 periodEnv
      vecs =
        asteriskGaussianFull
          numR2Freqs
          thetaFreq
          rFreq
          alpha
          periodR2
          periodEnv
          stdTheta
          stdR
  gaussian2DF <-
    dftExecute plan planID . VU.convert . toUnboxed . computeS . makeFilter2D $
    gaussian2D
  rF <-
    dftExecute plan planID . VU.convert . toUnboxed . computeS . makeFilter2D $
    r
  let filterF = VS.zipWith (\x y -> x + y ^ 2) gaussian2DF rF
  asteriskGaussianF <- dftExecuteBatchP plan planID vecs
  fmap VS.concat .
    dftExecuteBatchP plan inversePlanID .
    parMap rdeepseq (VS.zipWith (*) filterF) $
    asteriskGaussianF
