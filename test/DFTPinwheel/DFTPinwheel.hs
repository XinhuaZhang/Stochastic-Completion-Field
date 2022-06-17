{-# LANGUAGE FlexibleContexts #-}
module DFTPinwheel where

import           Control.Monad                  as M
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           DFT.Plan
import           Filter.Utils
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.BlockCudaMatrix
import           FourierMethod.FourierSeries2D
import           Image.IO
import           Pinwheel.Base
import           Pinwheel.FourierSeries2D
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random
import           Utils.SimpsonRule
import Utils.Array
import Numeric.GSL.Special.Bessel
import           Utils.Distribution
import Text.Printf

main = do
  args@(deviceIDsStr:numPointsStr:numR2FreqStr:deltaStr:deltaReconStr:periodR2Str:angularFreqStr:radialFreqStr:sigmaStr:periodEnvelopeStr:stdStr:translationStr:rotationAngleStr:scaleFactorStr:numBatchStr:_) <-
    getArgs
  let deviceIDs = read deviceIDsStr :: [Int]
      numPoints = read numPointsStr :: Int
      numR2Freqs = read numR2FreqStr :: Int
      delta = read deltaStr :: Double
      deltaRecon = read deltaReconStr :: Double
      periodR2 = read periodR2Str :: Double
      angularFreq = read angularFreqStr :: Int
      radialFreq = read radialFreqStr :: Int
      sigma = read sigmaStr :: Double
      periodEnvelope = read periodEnvelopeStr :: Double
      std = read stdStr :: Double
      translation = read translationStr :: Int
      rotationAngle = read rotationAngleStr :: Int
      scaleFactor = read scaleFactorStr :: Double
      numBatch = read numBatchStr :: Int
      folderPath = "output/test/DFTPinwheel"
  -- removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let center = div numPoints 2
      pinwheel =
        computeUnboxedS . R.fromFunction (Z :. numPoints :. numPoints) $ \(Z :. x :. y) ->
          fourierMellinInv
            sigma
            -- (periodEnvelope -- * sqrt 2
            -- )
            angularFreq
            radialFreq
            ( deltaRecon * fromIntegral (x - center)
            , deltaRecon * fromIntegral (y - center))
      -- pinwheelMat =
      --   CuMat (numPoints ^ 2) 1 . CuVecHost . VU.convert . toUnboxed $ pinwheel
      centerFreq = div numR2Freqs 2
      -- pinwheelFreq =
      --   computeUnboxedS . R.fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. x :. y) ->
      --     fourierMellin
      --       sigma
      --       angularFreq
      --       radialFreq
      --       ( delta * fromIntegral (x - centerFreq)
      --       , delta * fromIntegral (y - centerFreq))
      -- pinwheelFreqMat =
      --   CuMat 1 (numR2Freqs ^ 2) . CuVecHost . VU.convert . toUnboxed $ -- . computeS . weightedArray
      --   pinwheelFreq
  -- invHarmonics <-
  --   createInverseHarmonicMatriesGPU ptxs 1 numPoints numR2Freqs periodR2 delta
  -- pinwheelCoef <-
  --   computeFourierCoefficientsR2Stream
  --     deviceIDs
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     delta
  --     1
  --     numBatch
  --     [pinwheelMat]
  -- pinwheelFreqSeries <-
  --   computeFourierSeriesR2
  --     deviceIDs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     invHarmonics
  --     [pinwheelFreqMat]
  plotImageRepaComplex
    (folderPath </> (printf "Pinwheel_%d_%d.png" angularFreq radialFreq)) .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . R.traverse pinwheel id $ \f idx@(Z :. i :. j) ->
    let x = (fromIntegral $ i - center)
        y = (fromIntegral $ j - center)
        r = sqrt $ x ^ 2 + y ^ 2
       -- if r <= 6
       --   then 0
       --   else
     in f idx
  -- plotImageRepaComplex (folderPath </> "PinwheelTest.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelTest
  -- plotImageRepaComplex (folderPath </> "PinwheelCoefficients.png") .
  --   ImageRepa 8 . computeS . R.traverse pinwheelCoef id $ \f idx@(Z :. _ :. i :. j) ->
  --   let x = fromIntegral $ i - div numR2Freqs 2
  --       y = fromIntegral $ j - div numR2Freqs 2
  --       r = sqrt $ x ^ 2 + y ^ 2
  --   in if r <= 0
  --        then 0
  --        else f idx
  -- plotImageRepaComplex (folderPath </> "PinwheelFrequency.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelFreq
  -- plotImageRepaComplex (folderPath </> "PinwheelFrequencySeries.png") .
  --   ImageRepa 8 $
  --   pinwheelFreqSeries
  initVec <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numPoints ^ 2) randomIO
  lock <- getFFTWLock
  plan <-
    fst <$> idft1dGPlan lock emptyPlan [numPoints, numPoints] [0, 1] initVec
  -- dftPinwheel <-
  --   dftExecute plan (DFTPlanID IDFT1DG [numPoints, numPoints] [0, 1]) .
  --   VU.convert . toUnboxed . computeUnboxedS . makeFilter2D $
  --   pinwheel
  -- plotImageRepaComplex (folderPath </> "DFTPinwheel.png") .
  --   ImageRepa 8 .
  --   computeUnboxedS .
  --   makeFilter2D .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) . VS.convert $
  --   dftPinwheel
  let centerR2Freq = div numR2Freqs 2
      gaussian2D =
        computeUnboxedS . fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i' :. j') ->
          let i = i' - centerR2Freq
              j = j' - centerR2Freq
              rho = (sqrt . fromIntegral $ i ^ 2 + j ^ 2) / periodR2 * 2 * pi
               -- rho ^ (abs angularFreq) *
           in (exp
                 (pi * fromIntegral (i ^ 2 + j ^ 2) /
                  ((-1) * periodR2 ^ 2 * std ^ 2)) /
               (2 * pi * std ^ 2)) :+
              0
      -- pinwheelFreqAnaticalNonGaussian =
      --   centerHollowArray numR2Freqs . computeUnboxedS $
      --   analyticalFourierCoefficients2
      --     numR2Freqs
      --     delta
      --     angularFreq
      --     radialFreq
      --     (-sigma) --(sigma - 1)
      --     periodR2
      --     (periodEnvelope * sqrt 2)
      pinwheelFreqAnatical' =
        R.map
          (\x ->
             x *
             exp
               (0 :+
                ((-1) * fromIntegral angularFreq * fromIntegral rotationAngle * 5 / 180 * pi))) .
        centerHollowArray numR2Freqs . computeUnboxedS $
        analyticalFourierCoefficients2
          numR2Freqs
          delta
          angularFreq
          radialFreq
          sigma --(sigma - 1)
          periodR2
          (periodEnvelope -- * sqrt 2
           )
      pinwheelFreqAnatical =
        R.traverse pinwheelFreqAnatical' id $ \f idx@(Z :. x' :. y') ->
          let center = div numR2Freqs 2
              x = fromIntegral $ x' - center
              y = fromIntegral $ y' - center
              r = sqrt $ x ^ 2 + y ^ 2
           in f idx * (exp (0 :+ 2 * pi / periodR2 * ((-1) * (x + y) * fromIntegral translation))) *
              ((scaleFactor :+ 0) **
               (sigma :+
                (-1) * (2 * pi) / (log periodEnvelope) *
                (fromIntegral radialFreq)))
      pinwheelFreqAnaticalMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freqs ^ 2) A.:. (1 :: Int)) . R.toList $
        pinwheelFreqAnatical
         -- pinwheelFreqAnaticalInner =
         --   centerHollowArray numR2Freqs . computeUnboxedS $
         --   analyticalFourierCoefficients2'
         --     numR2Freqs
         --     (1 / fromIntegral (div numR2Freqs 2))
         --     angularFreq
         --     radialFreq
         --     (-sigma)
         --     periodR2
         --     (periodEnvelope * sqrt 2)
         -- pinwheelFreqAnaticalInnerMat =
         --   A.use .
         --   A.fromList (A.Z A.:. (numR2Freqs ^ 2) A.:. (1 :: Int)) . R.toList $
         --   pinwheelFreqAnaticalInner
         -- pinwheelFreqSeriesAnatical =
         --   computeUnboxedS $
         --   analyticalFourierSeries1
         --     numPoints
         --     1 -- delta
         --     angularFreq
         --     radialFreq
         --     (sigma - 1)
         --     periodR2
  if numPoints < numR2Freqs
    then let diff = div (numR2Freqs - numPoints) 2
          in plotImageRepaComplex
               (folderPath </>
                (printf "FourierPinwheel_%d_%d.png" angularFreq radialFreq)) .
             ImageRepa 8 .
             computeS .
             extend (Z :. (1 :: Int) :. All :. All) .
             R.backpermute
               (Z :. numPoints :. numPoints)
               (\(Z :. i :. j) -> (Z :. i + diff :. j + diff)) $
             pinwheelFreqAnatical
    else plotImageRepaComplex
           (folderPath </>
            (printf "FourierPinwheel_%d_%d.png" angularFreq radialFreq)) .
         ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
         pinwheelFreqAnatical
  -- plotImageRepaComplex
  --   (folderPath </> "PinwheelCoefficientsAnalyticalNonGaussian.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelFreqAnaticalNonGaussian
  -- plotImageRepaComplex
  --   (folderPath </> "PinwheelCoefficientsAnalyticalInner.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelFreqAnaticalInner
  -- plotImageRepaComplex (folderPath </> "PinwheelFrequencySeriesAnalytical.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelFreqSeriesAnatical
  -- pinwheelFreqSeriesStream <-
  --   computeFourierSeriesR2Stream
  --     deviceIDs
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     delta
  --     numBatch
  --     [transposeCuMat pinwheelFreqMat]
  -- plotImageRepaComplex (folderPath </> "PinwheelFrequencySeriesStream.png") .
  --   ImageRepa 8 . computeS . R.traverse pinwheelFreqSeriesStream id $ \f idx@(Z :. _ :. i :. j) ->
  --   let x = fromIntegral $ i - div numPoints 2
  --       y = fromIntegral $ j - div numPoints 2
  --       r = sqrt $ x ^ 2 + y ^ 2
  --   in if r <= 0
  --        then 0
  --        else f idx
  pinwheelCoefficientsAnalyticalSeriesStream <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freqs
      numPoints
      1
      periodR2
      deltaRecon
      numBatch
      pinwheelFreqAnaticalMat
      -- [ transposeCuMat . createCuMat $
      --   [VU.convert . toUnboxed $ pinwheelFreqAnatical]
      -- ]
  plotImageRepaComplex
    (folderPath </>
     (printf
        "Pinwheel_FourierSeries_%d_%d_%d_%03d_%02d_%.1f.png"
        angularFreq
        radialFreq
        numR2Freqs
        translation
        rotationAngle
        scaleFactor)) .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    R.traverse (sumS . rotate3D $ pinwheelCoefficientsAnalyticalSeriesStream) id $ \f idx@(Z :. i :. j) ->
    let x = fromIntegral $ i - center
        y = fromIntegral $ j - center
        r = sqrt $ x ^ 2 + y ^ 2
       -- if r <= 6
       --   then 0
       --   else
     in f idx
  -- pinwheelCoefficientsAnalyticalIDFT <-
  --   dftExecute plan (DFTPlanID IDFT1DG [numPoints, numPoints] [0, 1]) .
  --   VU.convert . toUnboxed . computeS . makeFilter2D $
  --   pinwheelFreqAnatical
  -- plotImageRepaComplex (folderPath </> "PinwheelCozefficientsIDFT.png") .
  --   ImageRepa 8 .
  --   computeS .
  --   makeFilter2DInverse .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) . VS.convert $
  --   pinwheelCoefficientsAnalyticalIDFT
  -- pinwheelCoefficientsAnalyticalInnerSeriesStream <-
  --   computeFourierSeriesR2StreamAcc'
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     1
  --     periodR2
  --     1 -- (1 / fromIntegral (div numR2Freqs 2))
  --     numBatch
  --     pinwheelFreqAnaticalMat
  -- plotImageRepaComplex
  --   (folderPath </> "PinwheelCoefficientsAnalyticalInnverSeries.png") .
  --   ImageRepa 8 .
  --   computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
  --   pinwheelCoefficientsAnalyticalInnerSeriesStream
  -- plotImageRepaComplex
  --   (folderPath </> "PinwheelCoefficientsAnalyticalSumSeries.png") .
  --   ImageRepa 8 .
  --   computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
  --   R.zipWith
  --     (+)
  --     pinwheelCoefficientsAnalyticalSeriesStream
  --     pinwheelCoefficientsAnalyticalInnerSeriesStream
  --   computeS . R.traverse pinwheelCoefficientsAnalyticalSeriesStream id $ \f idx@(Z :. _ :. i :. j) ->
  --   let x = fromIntegral $ i - div numPoints 2
  --       y = fromIntegral $ j - div numPoints 2
  --       r = sqrt $ x ^ 2 + y ^ 2
  --   in if r <= 0
  --        then 0
  --        else f idx
  -- -- pinwheelCoefficientsSeriesStream <-
  --   computeFourierSeriesR2Stream
  --     deviceIDs
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     delta
  --     numBatch
  --     [transposeCuMat . createCuMat $ [VU.convert . toUnboxed $ pinwheelCoef]]
  -- plotImageRepaComplex (folderPath </> "PinwheelCoefficientsSeriesStream.png") .
  --   ImageRepa 8 . computeS . R.traverse pinwheelCoefficientsSeriesStream id $ \f idx@(Z :. _ :. i :. j) ->
  --   let x = fromIntegral $ i - div numPoints 2
  --       y = fromIntegral $ j - div numPoints 2
  --       r = sqrt $ x ^ 2 + y ^ 2
  --   in if r <= 0
  --        then 0
  --        else f idx
  -- let besselFilter =
  --       computeUnboxedS . fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
  --         let xFreq = fromIntegral $ i - centerFreq
  --             yFreq = fromIntegral $ j - centerFreq
  --             rho = sqrt $ xFreq ^ 2 + yFreq ^ 2
  --         in if rho == 0
  --              then 1
  --              else (bessel_J1 (2 * pi * rho / periodR2)) / rho :+ 0
  -- plotImageRepaComplex (folderPath </> "BesselFilter.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   besselFilter
  -- initVec1 <-
  --   (VS.fromList . L.map (\x -> x :+ 0)) <$>
  --   M.replicateM (numR2Freqs ^ 2) randomIO
  -- lock <- getFFTWLock
  -- plan <-
  --   fst <$>
  --   (dft1dGPlan lock emptyPlan [numR2Freqs, numR2Freqs] [0, 1] initVec1 >>= \(plan, vec) ->
  --      idft1dGPlan lock plan [numR2Freqs, numR2Freqs] [0, 1] vec)
  -- let besselVec =
  --       VU.convert . toUnboxed . computeS . makeFilter2D $ besselFilter
  --     pinwheelFreqAnaticalVec = VU.convert . toUnboxed $ pinwheelFreqAnatical
  -- besselVecF <-
  --   dftExecute plan (DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]) besselVec
  -- pinwheelFreqAnaticalVecF <-
  --   dftExecute
  --     plan
  --     (DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1])
  --     pinwheelFreqAnaticalVec
  -- pinwheelFreqAnatical1 <-
  --   fmap (fromUnboxed (Z :. numR2Freqs :. numR2Freqs) . VS.convert) .
  --   dftExecute plan (DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1]) $
  --   VS.zipWith (*) besselVecF pinwheelFreqAnaticalVecF
  -- plotImageRepaComplex
  --   (folderPath </> "BesselFilterConveledPinwheelFreqAnatical.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelFreqAnatical1
  -- let pinwheelFreqAnatical2 =
  --       R.zipWith (-) pinwheelFreqAnatical pinwheelFreqAnatical1
  --     pinwheelFreqAnaticalMat2 =
  --       A.use .
  --       A.fromList (A.Z A.:. (numR2Freqs ^ 2) A.:. (1 :: Int)) . R.toList $
  --       pinwheelFreqAnatical2
  -- plotImageRepaComplex (folderPath </> "pinwheelFreqAnatical2.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelFreqAnatical2
  -- pinwheelCoefficientsAnalyticalSeriesStream2 <-
  --   computeFourierSeriesR2StreamAcc
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     1
  --     periodR2
  --     delta
  --     numBatch
  --     pinwheelFreqAnaticalMat2
  -- plotImageRepaComplex
  --   (folderPath </> "PinwheelCoefficientsAnalyticalSeries2.png") .
  --   ImageRepa 8 .
  --   computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
  --   pinwheelCoefficientsAnalyticalSeriesStream2
  -- dftPinwheel <-
  --   fmap (fromUnboxed (Z :. numR2Freqs :. numR2Freqs) . VS.convert) .
  --   dftExecute plan (DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]) .
  --   VU.convert . toUnboxed . computeS . makeFilter2D $
  --   pinwheel
  -- let zeroArr =
  --       fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
  --         if (sqrt . fromIntegral $ (i - centerFreq) ^ 2 + (j - centerFreq) ^ 2) <
  --            1
  --           then 0
  --           else 1
  -- plotImageRepaComplex (folderPath </> "DFTPinwheel.png") .
  --   ImageRepa 8 .
  --   computeS .
  --   extend (Z :. (1 :: Int) :. All :. All) .
  --   R.zipWith (*) zeroArr . makeFilter2D $
  --   dftPinwheel
  -- let normFunc arr =
  --       let m = VU.maximum . toUnboxed . computeS . R.map magnitude $ arr
  --       in R.map (/ (m :+ 0)) arr
  -- plotImageRepaComplex (folderPath </> "Diff.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   R.zipWith
  --     (\a b ->
  --        if magnitude b == 0
  --          then 0
  --          else a / b)
  --     (normFunc . R.zipWith (*) zeroArr . makeFilter2D $ dftPinwheel)
  --     (normFunc . R.zipWith (*) zeroArr $ pinwheelFreqAnatical)
