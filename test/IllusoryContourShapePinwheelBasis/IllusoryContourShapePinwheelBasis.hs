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
import           Utils.Time

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:tauStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdR2Str:stdStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:writeFlagStr:numIterationStr:shape2DStr:numThreadStr:_) <-
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
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      stdR2 = read stdR2Str :: Double
      std = read stdStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      numBatchOri = read numBatchOriStr :: Int
      batchSize = read batchSizeStr :: Int
      s = read sStr :: Double
      writeFlag = read writeFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      shape2D@(Points _ minDist _) = read shape2DStr :: Points Shape2D
  -- removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  hist <-
    if False
      then do
        printCurrentTime "read coefficients from file"
        decodeFile histFilePath
      else do
        printCurrentTime "Start computing coefficients..."
        sampleCartesian
          histFilePath
          folderPath
          ptxs
          numPoints
          delta
          numOrientation
          initScale
          thetaSigma
          tau
          threshold
          s
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
  printCurrentTime "Done"
  printCurrentTime "Start Convloution.."
  let coefficients =
        normalizeFreqArr' std phiFreqs rhoFreqs . getNormalizedHistogramArr $
        hist
      harmonicsArray =
        pinwheelFourierCoefficientsAnatical
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          (-s)
          periodR2
      points =
        L.map (\(x, y) -> Point (fromIntegral x) (fromIntegral y) 0 1) . getShape2DIndexList . makeShape2D $
        shape2D
      initialSourceDistribution =
        computeInitialDistributionPowerMethodPinwheelBasis'
          numR2Freq
          periodR2
          phiFreq
          rhoFreq
          points
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPointsRecon
      numPointsRecon
      numR2Freq
      (2 * thetaFreq + 1)
      (2 * scaleFreq + 1)
  dftBias <- computeBiasPinwheelBasis plan numR2Freq periodR2 points
  let biasMat =
        A.use . A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VU.toList $
        computeBiasPinwheelBasis1 plan numR2Freq periodR2 points
  print points
  biasR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      1
      biasMat
  plotImageRepa (folderPath </> "Bias.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt .
    toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    biasR2
  computeContourPinwheelBasis
    plan
    folderPath
    ptxs
    writeFlag
    coefficients
    harmonicsArray
    dftBias
    numIteration
    numBatchR2
    numPointsRecon
    numR2Freq
    deltaRecon
    periodR2
    (2 * thetaFreq + 1)
    (2 * scaleFreq + 1)
    initialSourceDistribution
