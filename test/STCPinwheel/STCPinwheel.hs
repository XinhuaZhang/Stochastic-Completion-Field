{-# LANGUAGE BangPatterns #-}
module STCPinwheel where

import           Control.Monad                  as M
import           Control.Monad.Parallel         as MP
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.IArray              as IA
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           Filter.Utils
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.GreensFunction
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.BlockMatrixAcc
import           FourierMethod.FourierSeries2D
import           FourierPinwheel
import           Image.IO
import           Pinwheel.FourierSeries2D
import           STC                            hiding (convolve)
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
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
      radius = read radiusStr :: Double
      periodEnvelope = periodR2 * sqrt 2
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  printCurrentTime "Start computing coefficients..."
  hist <-
    if flag
      then decodeFile histFilePath
      else sampleCartesian
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
  let coefficients = getNormalizedHistogramArr $ hist
  -- coefficients <- hollowCoefficients plan radius periodEnvelope coefficients'
      -- coefficients = getNormalizedHistogramArr $ hist
        -- normalizeFreqArr' std phiFreqs rhoFreqs .
  let !harmonicsArray =
        createHarmonics
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          (-s)
          periodR2
          (periodR2 * sqrt 2)
          coefficients
      harmonicsArray1 =
        createHarmonics
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          (-1.25)
          periodR2
          (periodR2 * sqrt 2)
          coefficients
      !initialSourceDistribution =
        computeInitialDistributionFourierPinwheel
          numR2Freq
          periodR2
          (periodR2 * sqrt 2)
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          initSource
      !initialSinkDistribution =
        computeInitialDistributionFourierPinwheel
          numR2Freq
          periodR2
          (periodR2 * sqrt 2)
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          initSink
      -- harmonicsArray' =
      --   centerHollow numR2Freq $
      --   pinwheelFourierCoefficientsAnatical
      --     numR2Freq
      --     phiFreq
      --     rhoFreq
      --     thetaFreq
      --     scaleFreq
      --     (-s)
      --     periodR2
      --     periodEnvelope
      -- harmonicsArray1' =
      --   -- L.map (centerHollowVector numR2Freq) $
      --   pinwheelFourierCoefficientsAnaticalList1
      --     numR2Freq
      --     phiFreq
      --     rhoFreq
      --     thetaFreq
      --     scaleFreq
      --     (-s)
      --     periodR2
      --     periodEnvelope
      -- harmonicsArray2' =
      --   -- L.map (centerHollowVector numR2Freq) $
      --   pinwheelFourierCoefficientsAnaticalList2
      --     numR2Freq
      --     phiFreq
      --     rhoFreq
      --     thetaFreq
      --     scaleFreq
      --     (-s)
      --     periodR2
      --     periodEnvelope
      -- initialSourceDistribution =
      --   computeInitialDistributionFull'
      --     numR2Freq
      --     periodR2
      --     phiFreq
      --     rhoFreq
      --     initSource
      -- initialSourceDistribution =
      --   computeInitialDistributionPowerMethodPinwheelBasis'
      --     numR2Freq
      --     (-s)
      --     periodR2
      --     phiFreq
      --     rhoFreq
      --     initSource
      -- initialSinkDistribution =
      --   computeInitialDistributionFull'
      --     numR2Freq
      --     periodR2
      --     phiFreq
      --     rhoFreq
      --     initSink
  -- source <-
  --   convoluvePinhweelBasisTest'
  --     coefficients
  --     harmonicsArray1'
  --     harmonicsArray2'
  --     initialSourceDistribution
  -- sink <-
  --   timeReversal' <$>
  --   convoluvePinhweelBasisTest'
  --     coefficients
  --     harmonicsArray1'
  --     harmonicsArray2'
  --     initialSinkDistribution
  -- let source =
  --       convoluvePinhweelBasis'
  --         coefficients
  --         harmonicsArray'
  --         initialSourceDistribution
  --     sink =
  --       timeReversal' $
  --       convoluvePinhweelBasis'
  --         coefficients
  --         harmonicsArray'
  --         initialSinkDistribution
  --     sourceMat =
  --       A.transpose . A.use . fromDFTArray . getDFTArrayVector $ source
  --     sinkMat = A.transpose . A.use . fromDFTArray . getDFTArrayVector $ sink-
  source <- convolve harmonicsArray initialSourceDistribution
  source1 <- convolve harmonicsArray1 initialSourceDistribution
  let numLogPolarFreq = ((2 * scaleFreq + 1) * (2 * thetaFreq + 1))
      sourceMat = A.transpose . toMatrixAcc $ source
      sourceMat1 = A.transpose . toMatrixAcc $ source1
      sourceMat2 =
        A.transpose . toMatrixAcc $
        parZipWithFPArray (VS.zipWith (-)) source source1
  sourceR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      numLogPolarFreq
      periodR2
      deltaRecon
      numBatchR2
      sourceMat
  plotImageRepa (folderPath </> "SourcePower.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    sourceR2
  source1R2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      numLogPolarFreq
      periodR2
      deltaRecon
      numBatchR2
      sourceMat1
  plotImageRepa (folderPath </> "SourcePower1.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    source1R2
  source2R2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      numLogPolarFreq
      periodR2
      deltaRecon
      numBatchR2
      sourceMat2
  plotImageRepa (folderPath </> "SourcePower2.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    source2R2
  sink <- convolve harmonicsArray initialSinkDistribution
  let sinkMat = A.transpose . toMatrixAcc $ sink
  -- completion <- completionField' plan source (timeReversal' sink)
  -- let completionMat = createCuMat . getDFTArrayVector $ completion
  -- print completionMat
  sinkR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      numLogPolarFreq
      periodR2
      deltaRecon
      numBatchR2
      sinkMat
  plotImageRepa (folderPath </> "SinkPower.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    sinkR2
  -- completionR2 <-
  --   computeFourierSeriesR2Stream
  --     deviceIDs
  --     ptxs
  --     numR2Freq
  --     numPointsRecon
  --     periodR2
  --     deltaRecon
  --     numBatchR2
  --     [transposeCuMat completionMat]
  printCurrentTime "Start ploting.."
  -- plotDFTArrayPower
  --   (folderPath </> "CoefficientsPower.png")
  --   numR2Freq
  --   numR2Freq
  --   source
  -- plotImageRepaComplex (folderPath </> "Coefficients.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) .
  --   VU.convert . L.foldl1' (VS.zipWith (+)) . getDFTArrayVector $
  --   source
  -- plotImageRepaComplex (folderPath </> "Source.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
  --   toUnboxed . sumS . rotate3D $
  --   sourceR2
  completionR2 <-
    completionFieldRepa
      plan
      (R.reshape
         (Z :. (L.length scaleFreqs) :. (L.length thetaFreqs) :. numPointsRecon :.
          numPointsRecon) $
       sourceR2)
      (timeReversalRepa thetaFreqs .
       R.reshape
         (Z :. (L.length scaleFreqs) :. (L.length thetaFreqs) :. numPointsRecon :.
          numPointsRecon) $
       sinkR2)
  plotImageRepa (folderPath </> "CompletionPower.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    completionR2
  printCurrentTime "Done"
  -- sparse
  let harmonicsArraySparse =
        computeHarmonicsArraySparse
          numPointsRecon
          deltaRecon
          numPointsRecon
          deltaRecon
          phiFreqs
          rhoFreqs
          thetaFreqs
          scaleFreqs
          -- (2 * pi)
          (log (periodR2 / 4))
          64
          s
      init =
        computeInitialDistribution'
          numPointsRecon
          numPointsRecon
          phiFreqs
          rhoFreqs
          0
          initSource
  harmonicsArrayDFT <-
    fmap (listArray (bounds harmonicsArraySparse)) .
    dftExecuteBatchP
      plan
      (DFTPlanID DFT1DG [numPointsRecon, numPointsRecon] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    harmonicsArraySparse
  sourceSparse <- convolve' Source plan coefficients harmonicsArrayDFT init
  plotDFTArrayPower
    (folderPath </> "SourcePower_Sparse.png")
    numPointsRecon
    numPointsRecon
    sourceSparse
