{-# LANGUAGE Strict #-}
module STC.Bias where

import           Array.UnboxedArray       as AU
import           Data.Array.Repa          as R
import           Data.Complex
import           Data.List                as L
import           Data.Vector.Storable     as VS
import           Data.Vector.Unboxed      as VU
import Data.Vector.Generic as VG
import           DFT.Plan
import           Filter.Utils
import           Pinwheel.FourierSeries2D
import           STC.Point
import           STC.Utils
import           Utils.List
import FourierPinwheel.AsteriskGaussian
import Filter.Utils
-- import           Data.Array.Accelerate.LLVM.PTX as A
-- import qualified Data.Array.Accelerate          as A
-- import FourierMethod.FourierSeries2D
import Image.IO
import Utils.Array
import System.FilePath
import Utils.Parallel
import FourierPinwheel.GaussianEnvelopePinwheel

{-# INLINE computeBias #-}
computeBias ::
     (Num a, Unbox a, VG.Vector vector a) => Int -> Int -> [Point] -> vector a
computeBias rows cols =
  let (minR, maxR) = computeRange rows
      (minC, maxC) = computeRange cols
   in VG.convert .
      toUnboxedVector .
      AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
      L.map (\(Point x y theta scale) -> ((round x, round y), 1))

computeBiasPinwheelBasis ::
     DFTPlan -> Int -> Double -> Double -> [Point] -> IO (VS.Vector (Complex Double))
computeBiasPinwheelBasis plan numR2Freq sigma period points =
  let r2Freqs = getListFromNumber' numR2Freq
      envelope =
        computeUnboxedS $
        analyticalFourierCoefficients2 numR2Freq 1 0 0 sigma period (period * sqrt 2)
  in dftExecute plan (DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]) .
     VS.convert .
     toUnboxed .
     computeS .
     makeFilter2D .
     R.zipWith (*) envelope .
     fromListUnboxed (Z :. numR2Freq :. numR2Freq) $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasPinwheelBasis1 ::
     DFTPlan -> Int -> Double -> Double -> [Point] -> VU.Vector (Complex Double)
computeBiasPinwheelBasis1 plan numR2Freq sigma period points =
  let r2Freqs = getListFromNumber' numR2Freq
      envelope =
        toUnboxed . computeS $
        analyticalFourierCoefficients2 numR2Freq 1 0 0 sigma period (period * sqrt 2)
  in VU.zipWith (*) envelope . VU.fromList $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasGaussian ::
     DFTPlan
  -> Int
  -> Double
  -> Double
  -> [Point]
  -> IO (VS.Vector (Complex Double))
computeBiasGaussian plan numR2Freq sigma period points =
  let r2Freqs = getListFromNumber' numR2Freq
      freqCenter = div numR2Freq 2
      envelope =
        fromFunction (Z :. numR2Freq :. numR2Freq) $ \(Z :. i :. j) ->
          (exp $
           (fromIntegral $ (i - freqCenter) ^ 2 + (j - freqCenter) ^ 2) *
           (sigma ^ 2) /
           (-2)) *
          (sigma ^ 2) /
          (2 * pi) :+
          0
  in dftExecute plan (DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]) .
     VS.convert .
     toUnboxed .
     computeS .
     makeFilter2D .
     R.zipWith (*) envelope . fromListUnboxed (Z :. numR2Freq :. numR2Freq) $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasGaussian1 ::
     DFTPlan -> Int -> Double -> Double -> [Point] -> VU.Vector (Complex Double)
computeBiasGaussian1 plan numR2Freq sigma period points =
  let r2Freqs = getListFromNumber' numR2Freq
      freqCenter = div numR2Freq 2
      envelope =
        toUnboxed . computeS . fromFunction (Z :. numR2Freq :. numR2Freq) $ \(Z :. i :. j) ->
          (exp $
           (fromIntegral $ (i - freqCenter) ^ 2 + (j - freqCenter) ^ 2) *
           (sigma ^ 2) /
           (-2)) *
          (sigma ^ 2) /
          (2 * pi) :+
          0
  in VU.zipWith (*) envelope . VU.fromList $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasFourierPinwheel ::
     DFTPlan
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Point]
  -> IO (VS.Vector (Complex Double), VS.Vector (Complex Double))
computeBiasFourierPinwheel plan numR2Freqs thetaFreq alpha periodR2 periodEnv radius stdTheta points = do
  asteriskGaussianVec <-
    asteriskGaussianLowPass
      plan
      numR2Freqs
      thetaFreq
      alpha
      periodR2
      periodEnv
      stdTheta
      radius
  let r2Freqs = getListFromNumber' numR2Freqs
      numThetaFreq = 2 * thetaFreq + 1
      shiftVec =
        VS.concat . L.replicate numThetaFreq . VS.fromList $
        [ L.foldl'
          (\b (Point x y theta scale) ->
             b + cis (-(freqX * x + freqY * y) * 2 * pi / periodR2))
          0
          points
        | freqY <- r2Freqs
        , freqX <- r2Freqs
        ]
      biasVec = VS.zipWith (*) asteriskGaussianVec shiftVec
      -- centerFreq = div numR2Freqs 2
      -- gaussian2D =
      --   VU.convert .
      --   toUnboxed . computeS . fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
      --     let i = i' - centerFreq
      --         j = j' - centerFreq
      --      in exp
      --           (pi * fromIntegral (i ^ 2 + j ^ 2) /
      --            ((-1) * periodR2 ^ 2 * stdR2 ^ 2)) /
      --         (2 * pi * stdR2 ^ 2) :+
      --         0
      -- zeroVec = VS.replicate (numR2Freqs ^ 2) 0
      -- shiftVec =
      --   VS.concat .
      --   L.map
      --     (\af ->
      --        if af == 0
      --          then VS.zipWith (*) gaussian2D . VS.fromList $
      --               [ L.foldl'
      --                 (\b (Point x y theta scale) ->
      --                    b + cis (-(freqX * x + freqY * y) * 2 * pi / periodR2))
      --                 0
      --                 points
      --               | freqY <- r2Freqs
      --               , freqX <- r2Freqs
      --               ]
      --          else zeroVec) $
      --   [-thetaFreq .. thetaFreq]
  biasVecF <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreq, numR2Freqs, numR2Freqs] [0, 1, 2]) .
    VU.convert .
    toUnboxed .
    computeS .
    makeFilter3D .
    fromUnboxed (Z :. numThetaFreq :. numR2Freqs :. numR2Freqs) . VS.convert $
    -- shiftVec
    biasVec
  return (biasVec, biasVecF)
  -- return (shiftVec, biasVecF)
  
computeBiasFourierPinwheelFull ::
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
  -> VU.Vector (Complex Double)
  -> [Point]
  -> IO (VS.Vector (Complex Double), VS.Vector (Complex Double))
computeBiasFourierPinwheelFull plan numR2Freqs thetaFreq rFreq alpha periodR2  radius stdTheta stdR stdR2 asteriskGaussianVec points = do
  let r2Freqs = getListFromNumber' numR2Freqs
      shiftArr =
        fromListUnboxed (Z :. numR2Freqs :. numR2Freqs) .
        parMap
          rdeepseq
          (\(freqY, freqX) ->
             L.foldl'
               (\b (Point x y _ _) ->
                  b + cis (-(freqX * x + freqY * y) * 2 * pi / periodR2))
               0
               points) $
        [(freqY, freqX) | freqY <- r2Freqs, freqX <- r2Freqs]
      asteriskGaussianArr =
        fromUnboxed
          (Z :. rFreq :. thetaFreq :. numR2Freqs :. numR2Freqs)
          asteriskGaussianVec
  biasArr <-
    computeUnboxedP . R.traverse2 asteriskGaussianArr shiftArr const $ \fA fS idx@(Z :. _ :. _ :. i :. j) ->
      fA idx * fS (Z :. i :. j)
  biasVecF <-
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [rFreq, thetaFreq, numR2Freqs, numR2Freqs]
         [0, 1, 2, 3]) .
    VU.convert . toUnboxed . computeS . makeFilter4D $
    biasArr
  return
    ( VU.convert . toUnboxed . computeUnboxedS . makeFilter2D $ biasArr
    , biasVecF)
    
computeBiasFourierPinwheelKoffkaCorss ::
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
  -> VU.Vector (Complex Double)
  -> [Point]
  -> IO (VS.Vector (Complex Double), VS.Vector (Complex Double))
computeBiasFourierPinwheelKoffkaCorss plan numR2Freqs thetaFreq rFreq alpha periodR2  radius stdTheta stdR stdR2 asteriskGaussianVec points = do
  let r2Freqs = getListFromNumber' numR2Freqs
      thetaFreqs = getListFromNumber' thetaFreq
      halfPI = 0.5 * pi
      centerFreqR2 = div numR2Freqs 2
      centerThetaFreq = div thetaFreq 2
      shiftArr =
        fromFunction (Z :. thetaFreq :. numR2Freqs :. numR2Freqs) $ \(Z :. tf :. yf :. xf) ->
          let freqTheta = fromIntegral (tf - centerThetaFreq)
              freqY = fromIntegral (yf - centerFreqR2)
              freqX = fromIntegral (xf - centerFreqR2)
           in L.foldl'
                (\b (Point x y theta scale) ->
                   let phi =
                         if abs x < abs y
                           then atan2 0 x
                           else if abs x > abs y
                                  then atan2 y 0
                                  else error
                                         "computeBiasFourierPinwheelKoffkaCorss: x == y"
                    in b +
                       cis (-(freqX * x + freqY * y) * 2 * pi / periodR2) *
                       (cis (-freqTheta * (phi + halfPI)) +
                        cis (-freqTheta * (phi - halfPI))))
                0
                points
      asteriskGaussianArr =
        fromUnboxed
          (Z :. rFreq :. thetaFreq :. numR2Freqs :. numR2Freqs)
          asteriskGaussianVec
  biasArr <-
    computeUnboxedP . R.traverse2 asteriskGaussianArr shiftArr const $ \fA fS idx@(Z :. _ :. k :. i :. j) ->
      fA idx * fS (Z :. k :. i :. j)
  biasVecF <-
    (computeUnboxedP . makeFilter4D $ biasArr) >>=
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [rFreq, thetaFreq, numR2Freqs, numR2Freqs]
         [0, 1, 2, 3]) .
    VU.convert . toUnboxed
  biasVec <- VU.convert . toUnboxed <$> computeUnboxedP (makeFilter2D biasArr)
  return (biasVec, biasVecF)

computeBiasFourierPinwheelDiscrete ::
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
  -> Double
  -> [Point]
  -> IO (VS.Vector (Complex Double), VS.Vector (Complex Double))
computeBiasFourierPinwheelDiscrete plan numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv radius stdTheta stdR stdR2 points = do
  let biasVec =
        VU.convert .
        toUnboxed .
        computeS . makeFilter2D . fromUnboxed (Z :. numR2Freqs :. numR2Freqs) $
        computeBias numR2Freqs numR2Freqs points
      bias = VS.concat . L.replicate (2 * thetaFreq + 1) $ biasVec
      zeroVec = VS.replicate (VS.length biasVec) 0
  biasVecF <-
    dftExecute plan (DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]) biasVec
  let biasF =
        VS.concat .
        L.map
          (\angularFreq ->
             if angularFreq == 0
               then biasVecF
               else zeroVec) $
        [-thetaFreq .. thetaFreq]
  return (bias, biasF)
  

computeBiasFourierPinwheelDiscrete' ::
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
  -> Double
  -> [Point]
  -> IO (VS.Vector (Complex Double), VS.Vector (Complex Double))
computeBiasFourierPinwheelDiscrete' plan numR2Freqs thetaFreq rFreq alpha periodR2 periodEnv radius stdTheta stdR stdR2 points = do
  let r2Freqs = getListFromNumber' numR2Freqs
      shiftArr =
        fromListUnboxed (Z :. numR2Freqs :. numR2Freqs) .
        parMap
          rdeepseq
          (\(freqY, freqX) ->
             L.foldl'
               (\b (Point x y theta scale) ->
                  b + cis (-(freqX * (fromIntegral (round x :: Int)) + freqY * (fromIntegral (round y :: Int))) * 2 * pi / periodR2))
               0
               points) $
        [(freqY, freqX) | freqY <- r2Freqs, freqX <- r2Freqs]
      shiftVec =
        VU.convert . toUnboxed . computeUnboxedS . makeFilter2D $ shiftArr
      centerFreq = div numR2Freqs 2
      envelope =
        fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
          if i == centerFreq && j == centerFreq
            then 0
            else let r =
                       sqrt . fromIntegral $
                       (i - centerFreq) ^ 2 + (j - centerFreq) ^ 2
                  in (1 / (r * periodR2)) :+ 0
      envelopeShiftVec =
        VU.convert . toUnboxed . computeS . makeFilter2D $
        shiftArr -- *^ (centerHollowArray numR2Freqs envelope)
      zeroVec = VS.replicate (VS.length shiftVec) 0
  biasVec <-
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1])
      envelopeShiftVec
  let bias = VS.concat . L.replicate (2 * thetaFreq + 1) $ biasVec
      biasF =
        VS.concat .
        L.map
          (\angularFreq ->
             if angularFreq == 0
               then shiftVec
               else zeroVec) $
        [-thetaFreq .. thetaFreq]
  return (bias, biasF)

-- computeBiasFourierPinwheel ::
--      DFTPlan
--   -> Int
--   -> Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> [Point]
--   -> FilePath
--   -> [PTX]
--   -> Int
--   -> Int
--   -> IO (VS.Vector (Complex Double))
-- computeBiasFourierPinwheel plan numR2Freqs thetaFreq alpha periodR2 periodEnv stdR2 stdTheta points folderPath ptxs numPoints numBatch = do
--   let r2Freqs = getListFromNumber' numR2Freqs
--       numThetaFreq = 2 * thetaFreq + 1
--       gaussianVec =
--         asteriskGaussian1 numR2Freqs thetaFreq alpha periodR2 periodEnv stdR2
--       shiftCoef =
--         [ L.foldl'
--           (\b (Point x y theta scale) ->
--              b + cis (-(freqX * x + freqY * y) * 2 * pi / periodR2))
--           0
--           points
--         | freqY <- r2Freqs
--         , freqX <- r2Freqs
--         ]
--       shiftVec = VS.concat . L.replicate numThetaFreq . VS.fromList $ shiftCoef
--       biasVec = VS.zipWith (*) (VS.concat gaussianVec) shiftVec
--       biasMat =
--         A.transpose .
--         A.use .
--         A.fromList (A.Z A.:. numThetaFreq A.:. (numR2Freqs ^ 2)) . VS.toList $
--         biasVec
--   biasR2 <-
--     computeFourierSeriesR2StreamAcc
--       ptxs
--       numR2Freqs
--       numPoints
--       numThetaFreq
--       periodR2
--       1
--       numBatch
--       biasMat
--   plotImageRepa (folderPath </> "AsteriskBias.png") .
--     ImageRepa 8 .
--     fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
--     VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
--     biasR2
--   let centerFreq = div numR2Freqs 2
--       shiftArr = fromListUnboxed (Z :. numR2Freqs :. numR2Freqs) shiftCoef
--       gaussian2D =
--         computeUnboxedS . fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
--           let i = i' - centerFreq
--               j = j' - centerFreq
--            in (exp
--                  ((pi * (fromIntegral $ i ^ 2 + j ^ 2)) /
--                   ((-1) * periodR2 ^ 2 * stdR2 ^ 2))) /
--               (2 * pi * stdR2 ^ 2) :+
--               0
--       gaussian2DMat =
--         A.use .
--         A.fromList (A.Z A.:. (numR2Freqs ^ 2) A.:. (1 :: Int)) .
--         R.toList . R.zipWith (*) shiftArr $
--         gaussian2D
--   gaussian2DR2 <-
--     computeFourierSeriesR2StreamAcc
--       ptxs
--       numR2Freqs
--       numPoints
--       1
--       periodR2
--       1
--       numBatch
--       gaussian2DMat
--   plotImageRepa (folderPath </> "Gaussian2D.png") .
--     ImageRepa 8 .
--     fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
--     VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
--     gaussian2DR2
--   let planID = DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]
--       inversePlanID = DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1]
--   gaussian2DF <-
--     dftExecute plan planID . VU.convert . toUnboxed . computeS . makeFilter2D $
--     gaussian2D
--   asteriskGaussianF <- dftExecuteBatchP plan planID gaussianVec
--   outputVec <-
--     fmap VS.concat .
--     dftExecuteBatchP plan inversePlanID .
--     parMap rdeepseq (VS.zipWith (*) gaussian2DF) $
--     asteriskGaussianF
--   let outputMat =
--         A.transpose .
--         A.use .
--         A.fromList (A.Z A.:. numThetaFreq A.:. (numR2Freqs ^ 2)) . VS.toList $
--         outputVec
--   outputR2 <-
--     computeFourierSeriesR2StreamAcc
--       ptxs
--       numR2Freqs
--       numPoints
--       numThetaFreq
--       periodR2
--       1
--       numBatch
--       outputMat
--   plotImageRepa (folderPath </> "AsteriskGaussian.png") .
--     ImageRepa 8 .
--     fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
--     VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
--     outputR2
--   dftExecute
--     plan
--     (DFTPlanID DFT1DG [numThetaFreq, numR2Freqs, numR2Freqs] [0, 1, 2]) .
--     VU.convert .
--     toUnboxed .
--     computeS .
--     makeFilter3D .
--     fromUnboxed (Z :. numThetaFreq :. numR2Freqs :. numR2Freqs) . VS.convert $
--     -- shiftVec
    -- VS.zipWith (*) outputVec shiftVec



-- zeroVec = VS.replicate (numR2Freqs ^ 2) 0
-- shiftVec =
--   VS.concat .
--   L.map
--     (\angularFreq ->
--        if angularFreq == 0
--          then VS.fromList $
--               [ L.foldl'
--                 (\b (Point x y _ _) ->
--                    b + cis (-(freqX * x + freqY * y) * 2 * pi / periodR2))
--                 0
--                 points
--               | freqY <- r2Freqs
--               , freqX <- r2Freqs
--               ]
--          else zeroVec) $
--   [-thetaFreq .. thetaFreq]
