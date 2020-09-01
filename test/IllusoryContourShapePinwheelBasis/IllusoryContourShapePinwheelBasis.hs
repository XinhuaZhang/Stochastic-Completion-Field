module IllusoryContourShapePinwheelBasis where

import           Control.Monad                  as M
import           Control.Monad.Parallel         as MP
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           FokkerPlanck
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
import           Utils.Time
import           Utils.List
import Data.Array.IArray as IA
import FourierPinwheel

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:tauStr:numTrailsStr:deltaTStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdR2Str:stdThetaStr:stdRStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:writeFlagStr:numIterationStr:shape2DStr:radiusStr:numThreadStr:_) <-
    getArgs
  let deviceIDs = read deviceIDsStr :: [Int]
      numPoints = read numPointsStr :: Int
      delta = read deltaStr :: Double
      threshold = read thresholdStr :: Double
      numPointsRecon = read numPointsReconStr :: Int
      deltaRecon = read deltaReconStr :: Double
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      tau = read tauStr :: Double
      numTrails = read numTrailsStr :: Int
      deltaT = read deltaTStr :: Double
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
      initScale = read initScaleStr :: Double
      initDist = read initDistStr :: [(Double, Double, Double, Double)]
      initPoints = L.map (\(x, y, t, s) -> Point x y t s) initDist
      initSource = [L.head initPoints]
      initSink = [L.last initPoints]
      numThread = read numThreadStr :: Int
      folderPath = "output/test/IllusoryContourShapePinwheelBasis"
      halfLogPeriod = log maxScale
      stdR2 = read stdR2Str :: Double
      stdTheta = read stdThetaStr :: Double
      stdR = read stdRStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      numBatchOri = read numBatchOriStr :: Int
      batchSize = read batchSizeStr :: Int
      s = read sStr :: Double
      writeFlag = read writeFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      shape2D@(Points _ minDist shape) = read shape2DStr :: Points Shape2D
      radius = read radiusStr :: Double
      periodEnv = periodR2 * sqrt 2 -- ^ 2 * 2
      maxScale =  sqrt periodEnv -- periodR2 / 2 * sqrt 2
  -- removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  hist <-
    if flag
      then do
        printCurrentTime "read coefficients from file"
        decodeFile histFilePath
      else do
        printCurrentTime "Start computing coefficients..."
        initialise []
        devs <- M.mapM device deviceIDs
        ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
        ptxs <- M.mapM createTargetFromContext ctxs
        runMonteCarloFourierCoefficientsGPU
          deviceIDs
          numThread
          numTrails
          batchSize
          thetaSigma
          scaleSigma
          0
          s
          tau
          deltaT
          phiFreqs
          rhoFreqs
          thetaFreqs
          scaleFreqs          
          periodEnv
          histFilePath
        -- sampleCartesian
        --   histFilePath
        --   folderPath
        --   ptxs
        --   numPoints
        --   periodEnv
        --   delta
        --   numOrientation
        --   initScale
        --   thetaSigma
        --   tau
        --   threshold
        --   s
        --   phiFreq
        --   rhoFreq
        --   thetaFreq
        --   scaleFreq
  printCurrentTime "Done"
  printCurrentTime "Start Convloution.."
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPointsRecon
      numPointsRecon
      numR2Freq
      (2 * thetaFreq + 1)
      (2 * scaleFreq + 1)
      (2 * phiFreq + 1)
      (2 * rhoFreq + 1)
  printCurrentTime "Compute DFT Plan done."
  let coefficients = getNormalizedHistogramArr hist
  -- coefficients <- hollowCoefficients plan radius periodEnv coefficients'
      -- coefficients
      --   -- normalizeFreqArr stdA stdR phiFreqs rhoFreqs .
      --  = getNormalizedHistogramArr $ hist
  harmonicsArray <-
    createHarmonics
      numR2Freq
      phiFreq
      rhoFreq
      thetaFreq
      scaleFreq
      (-s)
      periodR2
      periodEnv
      coefficients
  let points =
        L.map (\(x, y) -> Point x y 0 1) . getShape2DIndexList' . makeShape2D $
        shape2D
  (bias, dftBias) <- --Full
  -- (dftBias, bias) <- -- Discrete
    computeBiasFourierPinwheelFull
      plan
      numR2Freq
      thetaFreq
      scaleFreq
      (-s)
      periodR2
      periodEnv
      radius
      stdTheta
      stdR
      stdR2
      points
  let initDist =
        computeInitialDistributionPowerMethodFourierPinwheelFull
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          bias -- dftBias
  -- initDist <- multiplyRFunction plan periodEnv initDist'
  computeContourFourierPinwheel
    plan
    folderPath
    writeFlag
    harmonicsArray
    dftBias -- bias
    numIteration
    numBatchR2
    numPointsRecon
    deltaRecon
    periodR2
    initDist
    (show . circleRadius $ shape)
  -- let initMat = A.transpose . toMatrixAcc $ initDist
  -- initR2 <-
  --   computeFourierSeriesR2StreamAcc
  --     ptxs
  --     (getFPArrayNumXFreq initDist)
  --     numPoints
  --     (getFPArrayNumRFreq initDist * getFPArrayNumThetaFreq initDist)
  --     periodR2
  --     delta
  --     numBatchR2
  --     initMat
  -- plotImageRepa (folderPath </> "Init.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
  --   VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
  --   initR2
