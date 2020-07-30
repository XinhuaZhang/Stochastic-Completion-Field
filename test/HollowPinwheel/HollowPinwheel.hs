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
  -- let a = read aStr :: Double
  --     b = read bStr :: Double
  --     delta = read deltaStr :: Double
  --     radialFreq = read radialFreqStr :: Int
  --     x = envelopIntegral a b delta (-s) periodR2 radialFreq
  --     xs =
  --       parMap
  --         rdeepseq
  --         (\y -> magnitude $ envelopIntegral a y delta (-s) periodR2 radialFreq)
  --         [(a + 3) .. b]
  --     m = round $ (b - a) / delta
  --     idxs = [(a + fromIntegral i * delta) | i <- [1 .. m]]
  --     ys =
  --       parMap
  --         rdeepseq
  --         (\i -> magnitude $ envelopIntegral2D numR2Freq 1 (-s) i radialFreq)
  --         idxs
  -- -- printEnvelopIntegral a b delta (-s) periodR2 radialFreq
  -- printf
  --   "Integral of r^(%.1f + i%d) in (%.1f, %.1f)\n%s %f\n"
  --   (-s)
  --   radialFreq
  --   a
  --   b
  --   (show x)
  --   (magnitude x)
  -- toFile def (folderPath </> "sum.png") $ do
  --   layout_title .=
  --     printf
  --       "r^(%.1f + i 2pi %d / log(%.1f)) (%.1f, %.1f) delta = %f"
  --       (-s)
  --       radialFreq
  --       (periodR2 / 2)
  --       a
  --       b
  --       delta
  --   plot (line "" [L.zip [(a + 3) .. b] xs])
  -- toFile def (folderPath </> "sum2D.png") $ do
  --   layout_title .= ""
  --   plot (line "" [L.zip idxs ys])
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
      -- hpf = laplacianHighPassFilter numR2Freq
      -- hpf = idealHighPassFilter1 numR2Freq
      -- lpf = idealLowPassFilter radius numR2Freq
      hpf = gaussianLowPassFilter std numR2Freq
      -- lpf = laplacianLowPassFilter std numR2Freq
      -- idealHPF = idealHighPassFilter radius numR2Freq
      -- gaussianHPF = gaussianHighPassFilter1 radius std numR2Freq
      ss = VS.sum (harmonicsArray' IA.! idx)
      -- hpf =
      --   VS.replicate (numR2Freq ^ 2) (ss / (fromIntegral (numR2Freq ^ 2) :+ 0))
      hpfMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        hpf
      hpfR2 =
        R.fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
          if i == div numPoints 2 && i == j
            then ss
            else 0
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPoints
      numPoints
      numR2Freq
      (2 * thetaFreq + 1)
      (2 * scaleFreq + 1)
  -- harmonicsArray <- convolveFrequency1 plan numR2Freq hpf harmonicsArray'
  let harmonicsArray = centerHollow numR2Freq harmonicsArray'
  let weightArr =
        computeWeightArrFromListOfShape [numR2Freq, numR2Freq] :: R.Array D DIM2 (Complex Double)
      -- simpsonWeight = VU.convert . toUnboxed . computeS $ weightArr
      pinwheelMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        (harmonicsArray' IA.! idx)
      pinwheelHPMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        (harmonicsArray IA.! idx)
  printf "ss = %f\n" (magnitude ss)
  printf "mag = %f\n" (magnitude . VS.sum $ (harmonicsArray' IA.! idx))
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
  let centerVal =
        (pinwheelR2 R.!
         (Z :. (0 :: Int) :. (div numPoints 2) :. (div numPoints 2)))
  print centerVal
  print . magnitude $ centerVal
  plotImageRepaComplex (folderPath </> "Pinwheel_HighPass.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    pinwheelHPR2
  plotImageRepaComplex (folderPath </> "Pinwheel_Coef.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    (harmonicsArray' IA.! idx)
  printf
    "filter center: %f\n" 
    (magnitude $
     (hpfR2' R.! (Z :. (0 :: Int) :. (div numPoints 2) :. (div numPoints 2))))
  plotImageRepaComplex (folderPath </> "HighPass.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    hpfR2'
  plotImageRepaComplex (folderPath </> "PinwheelHollow.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    sumS . rotate3D . R.zipWith (\x y -> y - x) hpfR2' $
    pinwheelR2
