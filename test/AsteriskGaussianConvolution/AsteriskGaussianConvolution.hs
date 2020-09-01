module AsteriskGaussianConvolution where

import           Control.Monad                  as M
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.GreensFunction
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.FourierSeries2D
import           FourierPinwheel
import           Image.IO
import           STC                            hiding (convolve)
import           System.Directory
import           System.Environment
import           System.FilePath
import           Utils.Array
import           Utils.List
import           Utils.Parallel
import           Utils.Time

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:tauStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdR2Str:stdStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:radiusStr:numThreadStr:_) <-
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
      folderPath = "output/test/AsteriskGaussianConvolution"
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      stdR2 = read stdR2Str :: Double
      std = read stdStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      numBatchOri = read numBatchOriStr :: Int
      batchSize = read batchSizeStr :: Int
      s = read sStr :: Double
      radius = read radiusStr :: Double
      periodEnvelope = periodR2^2 / 2 
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  hist <-
    if flag
      then do
        printCurrentTime "Read coefficients..."
        decodeFile histFilePath
      else do
        printCurrentTime "Start computing coefficients..."
        sampleCartesian
          histFilePath
          folderPath
          ptxs
          numPoints
          periodEnvelope
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
  let periodEnv = periodR2 * sqrt 2
      coefficients = getNormalizedHistogramArr hist
      harmonicsArray =
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
      asteriskGaussianVecs =
        asteriskGaussianFull
          numR2Freq
          thetaFreq
          scaleFreq
          (-s)
          periodR2
          periodEnv
          10
          10
      r2Freqs = getListFromNumber numR2Freq
      shiftFilter =
        fromListUnboxed
          (Z :. numR2Freq :. numR2Freq)
          [ L.foldl'
            (\b (Point x y theta scale) ->
               b +
               cis
                 (-(fromIntegral freqX * x + fromIntegral freqY * y) * 2 * pi /
                   periodR2))
            0
            initSource
          | freqY <- r2Freqs
          , freqX <- r2Freqs
          ]
      dftID = DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]
      idftID = DFTPlanID IDFT1DG [numR2Freq, numR2Freq] [0, 1]
  shiftFilterF <-
    dftExecute plan dftID . VU.convert . toUnboxed . computeS . makeFilter2D $
    shiftFilter
  asteriskGaussianF <- dftExecuteBatchP plan dftID asteriskGaussianVecs
  filteredAsteriskGaussian <-
    fmap VS.concat .
    dftExecuteBatchP plan idftID . parMap rdeepseq (VS.zipWith (*) shiftFilterF) $
    asteriskGaussianF
  let numThetaFreq = 2 * thetaFreq + 1
      numRFreq = 2 * scaleFreq + 1
      sourceDist =
        FPArray
          numR2Freq
          numR2Freq
          (2 * scaleFreq + 1)
          (2 * thetaFreq + 1)
          (2 * rhoFreq + 1)
          (2 * phiFreq + 1) .
        parMap
          rdeepseq
          (\radialFreq ->
             VU.convert .
             toUnboxed .
             computeS .
             R.slice
               (fromUnboxed
                  (Z :. numRFreq :. numThetaFreq :. numR2Freq :. numR2Freq) .
                VS.convert $
                filteredAsteriskGaussian) $
             (Z :. radialFreq :. All :. All :. All)) $
        [0 .. 2 * scaleFreq]
  source <- convolve harmonicsArray sourceDist
  let sourceMat = A.transpose . toMatrixAcc $ source
  sourceR2' <-
    computeFourierSeriesR2StreamAcc
      ptxs
      (getFPArrayNumXFreq source)
      numPointsRecon
      (getFPArrayNumRFreq source * getFPArrayNumThetaFreq source)
      periodR2
      delta
      numBatchR2
      sourceMat
  let sourceR2 =
        R.reshape
          (Z :. numRFreq :. numThetaFreq :. numPointsRecon :. numPointsRecon)
          sourceR2'
  plotImageRepa (folderPath </> "Source.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    sourceR2
