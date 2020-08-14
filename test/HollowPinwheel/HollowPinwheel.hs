module HollowPinwheel where

import           Control.Monad                  as M
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.IArray              as IA
import           Data.Array.Repa                as R
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.GreensFunction
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.BlockMatrixAcc
import           FourierMethod.FourierSeries2D
import           Image.IO
import           Pinwheel.FourierSeries2D
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.List
import           Utils.Time
import Graphics.Rendering.Chart.Easy
import Utils.Parallel
import Graphics.Rendering.Chart.Backend.Cairo
import Utils.SimpsonRule

main = do
  args@(deviceIDsStr:numPointsStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:stdStr:radiusStr:sStr:idxStr:numBatchStr:aStr:bStr:deltaStr:radialFreqStr:numThreadStr:_) <-
    getArgs
  let deviceIDs = read deviceIDsStr :: [Int]
      numPoints = read numPointsStr :: Int
      numR2Freq = read numR2FreqStr :: Int
      periodR2 = read periodR2Str :: Double
      phiFreq = read phiFreqsStr :: Int
      phiFreqs = L.map fromIntegral [-phiFreq .. phiFreq]
      rhoFreq = read rhoFreqsStr :: Int
      rhoFreqs = L.map fromIntegral [-rhoFreq .. rhoFreq]
      thetaFreq = read thetaFreqsStr :: Int
      thetaFreqs = L.map fromIntegral [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Int
      scaleFreqs = L.map fromIntegral [-scaleFreq .. scaleFreq]
      std = read stdStr :: Double
      radius = read radiusStr :: Double
      s = read sStr :: Double
      idx = read idxStr :: (Int, Int)
      numBatch = read numBatchStr :: Int
      numThread = read numThreadStr :: Int
      folderPath = "output/test/HollowPinwheel"
  createDirectoryIfMissing True folderPath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let harmonicsArray' =
        pinwheelFourierCoefficientsAnatical
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          (-s)
          periodR2
          (periodR2 * sqrt 2)
      hpf = idealLowPassFilter radius periodR2 numR2Freq
      hpfMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        hpf
  plotImageRepaComplex (folderPath </> "Filter_Coef.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    hpf
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPoints
      numPoints
      numR2Freq
      (2 * thetaFreq + 1)
      (2 * scaleFreq + 1)
  harmonicsArray <- convolveFrequency plan numR2Freq hpf harmonicsArray'
  -- let harmonicsArray = centerHollow numR2Freq harmonicsArray'
  let pinwheelMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        (harmonicsArray' IA.! idx)
      pinwheelHPMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        (harmonicsArray IA.! idx)
  plotImageRepaComplex (folderPath </> "FilteredPinhweel_Coef.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    (harmonicsArray IA.! idx)
  pinwheelR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPoints
      1
      periodR2
      1
      numBatch
      pinwheelMat
  pinwheelHPR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPoints
      1
      periodR2
      1
      numBatch
      pinwheelHPMat
  hpfR2' <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPoints
      1
      periodR2
      1
      numBatch
      hpfMat
  plotImageRepaComplex (folderPath </> "Pinwheel.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    pinwheelR2
  printf
    "pinwheel center: %f\n"
    (magnitude $
     (pinwheelR2 R.!
      (Z :. (0 :: Int) :. (div numPoints 2) :. (div numPoints 2 - 1))))
  plotImageRepaComplex (folderPath </> "Pinwheel_HighPass.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    pinwheelHPR2
  printf
    "lowpassed pinwheel center: %f\n"
    (magnitude $
     (pinwheelHPR2 R.!
      (Z :. (0 :: Int) :. (div numPoints 2) :. (div numPoints 2 - 1))))
  plotImageRepaComplex (folderPath </> "Pinwheel_Coef.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    (harmonicsArray' IA.! idx)
  printf
    "filter center: %f\n"
    (magnitude $
     (hpfR2' R.! (Z :. (0 :: Int) :. (div numPoints 2) :. (div numPoints 2 - 1))))
  plotImageRepaComplex (folderPath </> "HighPass.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    hpfR2'
  plotImageRepaComplex (folderPath </> "PinwheelHollow.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    sumS . rotate3D . R.zipWith (-) pinwheelR2 $
    pinwheelHPR2
