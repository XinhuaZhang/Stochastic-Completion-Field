{-# LANGUAGE Strict #-}
module DFTPinwheel where

import           Control.Monad                  as M
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

main = do
  args@(deviceIDsStr:numPointsStr:numR2FreqStr:deltaStr:periodR2Str:angularFreqStr:radialFreqStr:sigmaStr:_) <-
    getArgs
  let deviceIDs = read deviceIDsStr :: [Int]
      numPoints = read numPointsStr :: Int
      numR2Freqs = read numR2FreqStr :: Int
      delta = read deltaStr :: Double
      periodR2 = read periodR2Str :: Double
      angularFreq = read angularFreqStr :: Int
      radialFreq = read radialFreqStr :: Int
      sigma = read sigmaStr :: Double
      folderPath = "output/test/DFTPinwheel"
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let center = div numPoints 2
      pinwheel =
        computeUnboxedS . R.fromFunction (Z :. numPoints :. numPoints) $ \(Z :. x :. y) ->
          fourierMellin
            sigma
            angularFreq
            radialFreq
            (fromIntegral (x - center), fromIntegral (y - center))
      pinwheelMat =
        CuMat (numPoints ^ 2) 1 . CuVecHost . VU.convert . toUnboxed $ pinwheel
      centerFreq = div numR2Freqs 2
      pinwheelFreq =
        computeUnboxedS . R.fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. x :. y) ->
          fourierMellin
            sigma
            angularFreq
            radialFreq
            ( delta * fromIntegral (x - centerFreq)
            , delta * fromIntegral (y - centerFreq))
      pinwheelFreqMat =
        CuMat 1 (numR2Freqs ^ 2) . CuVecHost . VU.convert . toUnboxed $ -- . computeS . weightedArray
        pinwheelFreq
  -- invHarmonics <-
  --   createInverseHarmonicMatriesGPU ptxs 1 numPoints numR2Freqs periodR2 delta
  pinwheelCoef <-
    computeFourierCoefficientsR2
      deviceIDs
      ptxs
      numR2Freqs
      numPoints
      periodR2
      delta
      1
      1
      [pinwheelMat]
  -- pinwheelFreqSeries <-
  --   computeFourierSeriesR2
  --     deviceIDs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     invHarmonics
  --     [pinwheelFreqMat]
  -- plotImageRepaComplex (folderPath </> "Pinwheel.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheel
  -- plotImageRepaComplex (folderPath </> "PinwheelTest.png") .
  --   ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --   pinwheelTest
  plotImageRepaComplex (folderPath </> "PinwheelCoefficients.png") .
    ImageRepa 8 . computeS . R.traverse pinwheelCoef id $ \f idx@(Z :. _ :. i :. j) ->
    let x = fromIntegral $ i - div numR2Freqs 2
        y = fromIntegral $ j - div numR2Freqs 2
        r = sqrt $ x ^ 2 + y ^ 2
    in if r <= 2.5
         then 0
         else f idx
  plotImageRepaComplex (folderPath </> "PinwheelFrequency.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    pinwheelFreq
  -- plotImageRepaComplex (folderPath </> "PinwheelFrequencySeries.png") .
  --   ImageRepa 8 $
  --   pinwheelFreqSeries
  -- initVec <-
  --   (VS.fromList . L.map (\x -> x :+ 0)) <$>
  --   M.replicateM (numPoints ^ 2) randomIO
  -- lock <- getFFTWLock
  -- plan <-
  --   fst <$> dft1dGPlan lock emptyPlan [numPoints, numPoints] [0, 1] initVec
  -- dftPinwheel <-
  --   dftExecute plan (DFTPlanID DFT1DG [numPoints, numPoints] [0, 1]) .
  --   VU.convert . toUnboxed . computeUnboxedS . makeFilter2D $
  --   pinwheel
  -- plotImageRepaComplex (folderPath </> "DFTPinwheel.png") .
  --   ImageRepa 8 .
  --   computeUnboxedS .
  --   makeFilter2D .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) . VS.convert $
  --   dftPinwheel
  let pinwheelFreqAnatical =
        computeUnboxedS $
        analyticalFourierCoefficients1
          numR2Freqs
          1 -- delta
          angularFreq
          radialFreq
          (sigma - 1)
          periodR2
      pinwheelFreqSeriesAnatical =
        computeUnboxedS $
        analyticalFourierSeries1
          numPoints
          1 -- delta
          angularFreq
          radialFreq
          (sigma - 1)
          periodR2
  plotImageRepaComplex (folderPath </> "PinwheelCoefficientsAnalytical.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    pinwheelFreqAnatical
  plotImageRepaComplex (folderPath </> "PinwheelFrequencySeriesAnalytical.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    pinwheelFreqSeriesAnatical
  pinwheelFreqSeriesStream <-
    computeFourierSeriesR2Stream
      deviceIDs
      ptxs
      numR2Freqs
      numPoints
      periodR2
      delta
      2 -- numBatch
      [transposeCuMat pinwheelFreqMat]
  plotImageRepaComplex (folderPath </> "PinwheelFrequencySeriesStream.png") .
    ImageRepa 8 . computeS . R.traverse pinwheelFreqSeriesStream id $ \f idx@(Z :. _ :. i :. j) ->
    let x = fromIntegral $ i - div numPoints 2
        y = fromIntegral $ j - div numPoints 2
        r = sqrt $ x ^ 2 + y ^ 2
    in if r <= 2.5
         then 0
         else f idx
  pinwheelCoefficientsAnalyticalSeriesStream <-
    computeFourierSeriesR2Stream
      deviceIDs
      ptxs
      numR2Freqs
      numPoints
      periodR2
      delta
      2 -- numBatch
      [ transposeCuMat . createCuMat $
        [VU.convert . toUnboxed $ pinwheelFreqAnatical]
      ]
  plotImageRepaComplex
    (folderPath </> "PinwheelCoefficientsAnalyticalSeriesStream.png") .
    ImageRepa 8 .
    computeS . R.traverse pinwheelCoefficientsAnalyticalSeriesStream id $ \f idx@(Z :. _ :. i :. j) ->
    let x = fromIntegral $ i - div numPoints 2
        y = fromIntegral $ j - div numPoints 2
        r = sqrt $ x ^ 2 + y ^ 2
    in if r <= 2.5
         then 0
         else f idx
  pinwheelCoefficientsSeriesStream <-
    computeFourierSeriesR2Stream
      deviceIDs
      ptxs
      numR2Freqs
      numPoints
      periodR2
      delta
      2 -- numBatch
      [transposeCuMat . createCuMat $ [VU.convert . toUnboxed $ pinwheelCoef]]
  plotImageRepaComplex (folderPath </> "PinwheelCoefficientsSeriesStream.png") .
    ImageRepa 8 . computeS . R.traverse pinwheelCoefficientsSeriesStream id $ \f idx@(Z :. _ :. i :. j) ->
    let x = fromIntegral $ i - div numPoints 2
        y = fromIntegral $ j - div numPoints 2
        r = sqrt $ x ^ 2 + y ^ 2
    in if r <= 2.5
         then 0
         else f idx
