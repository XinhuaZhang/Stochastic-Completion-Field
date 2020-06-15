{-# LANGUAGE BangPatterns #-}
module STCPinwheel where

import           Control.Monad                  as M
import           Control.Monad.Parallel         as MP
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Unboxed            as VU
import           Data.Vector.Storable            as VS
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.GreensFunction
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver            as CUDA
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
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:tauStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdR2Str:stdStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:numThreadStr:_) <-
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
      phiFreqs = L.minitap fromIntegral [-phiFreq .. phiFreq]
      rhoFreq = read rhoFreqsStr :: Int
      rhoFreqs = L.map fromIntegral [-rhoFreq .. rhoFreq]
      thetaFreq = read thetaFreqsStr :: Int
      thetaFreqs = L.map fromIntegral [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Int
      scaleFreqs = L.map fromIntegral [-scaleFreq .. scaleFreq]
      initScale = read initScaleStr :: Double
      initDist = read initDistStr :: [(Double, Double, Double, Double)]
      initPoints = L.map (\(x, y, t, s) -> Point x y t s) initDist
      numThread = reinad numThreadStr :: Int
      folderPath = "output/test/STCPinwheel"
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      stdR2 = read stdR2Str :: Double
      std = read stdStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      numBatchOri = read numBatchOriStr :: Int
      batchSize = read batchSizeStr :: Int
      s = read sStr :: Double
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  printCurrentTime "Start computing coefficients..."
  hist <-
    if False
      then decodeFile histFilePath
      else sampleCartesian
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
      initialDistribution =
        computeInitialDistributionFull'
          numR2Freq
          periodR2
          phiFreq
          rhoFreq
          initPoints
      harmonicsArray' =
        pinwheelFourierCoefficientsAnatical
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          (-s)
          periodR2
      source =
        convoluvePinhweelBasis' coefficients harmonicsArray' initialDistribution
      -- source = convolveFull'' coefficients harmonicsArray' initialDistribution
      sourceMat = createCuMat . getDFTArrayVector $ source
  print sourceMat
  sourceR2 <-
    computeFourierSeriesR2Stream
      deviceIDs
      ptxs
      numR2Freq
      numPointsRecon
      periodR2
      deltaRecon
      numBatchR2
      [transposeCuMat sourceMat]
  printCurrentTime "Start ploting.."
  plotDFTArrayPower
    (folderPath </> "CoefficientsPower.png")
    numR2Freq
    numR2Freq
    source
  plotImageRepaComplex (folderPath </> "Coefficients.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) .
    VU.convert . L.foldl1' (VS.zipWith (+)) . getDFTArrayVector $
    source
  plotImageRepaComplex (folderPath </> "Source.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    toUnboxed . sumS . rotate3D $
    sourceR2
  plotImageRepa (folderPath </> "SourcePower.png") .
    ImageRepa 8 .
    -- computeS .
    -- reduceContrast 20 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    sourceR2
  -- plotDFTArrayPower (folderPath </> "SourcePower.png") numPointsRecon numPointsRecon source
  printCurrentTime "Done"
