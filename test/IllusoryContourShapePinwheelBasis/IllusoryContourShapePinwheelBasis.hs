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
import           Utils.List
import Data.Array.IArray as IA

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:tauStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdR2Str:stdAStr:stdRStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:writeFlagStr:numIterationStr:shape2DStr:radiusStr:numThreadStr:_) <-
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
      stdA = read stdAStr :: Double
      stdR = read stdRStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      numBatchOri = read numBatchOriStr :: Int
      batchSize = read batchSizeStr :: Int
      s = read sStr :: Double
      writeFlag = read writeFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      shape2D@(Points _ minDist _) = read shape2DStr :: Points Shape2D
      radius = read radiusStr :: Double
  -- removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  hist <-
    if flag
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
          (periodR2 * sqrt 2)
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
  let coefficients =
        normalizeFreqArr stdA stdR phiFreqs rhoFreqs . getNormalizedHistogramArr $
        hist
      harmonicsArray =
        centerHollow numR2Freq $
        pinwheelFourierCoefficientsAnatical
          numR2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          (-s)
          periodR2
          (periodR2 * sqrt 2)
      points =
        L.map (\(x, y) -> Point ( x) ( y) 0 1) .
        getShape2DIndexList' . makeShape2D $
        shape2D
      pinwheelMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        (harmonicsArray IA.! (0, 0))
  pinwheelR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      numBatchR2
      pinwheelMat
  plotImageRepaComplex (folderPath </> "Pinwheel.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    pinwheelR2
  dftBias <- computeBiasGaussian plan numR2Freq stdR2 periodR2 points
  let biasMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VU.toList $
        computeBiasGaussian1 plan numR2Freq stdR2 periodR2 points
  print points
  plotImageRepa (folderPath </> "BiasFreq.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VU.map magnitude $
    computeBiasGaussian1 plan numR2Freq stdR2 periodR2 points
  biasR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      numBatchR2
      biasMat
  plotImageRepa (folderPath </> "Bias.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    biasR2
  -- -- For debugging
  -- printCurrentTime "Discrete and sparse method"
  -- plan <-
  --   makePlan
  --     folderPath
  --     plan
  --     numPointsRecon
  --     numPointsRecon
  --     numPointsRecon
  --     (2 * thetaFreq + 1)
  --     (2 * scaleFreq + 1)
  -- let xs =
  --       L.map (\(a, b) -> (a, b)) . getShape2DIndexList' . makeShape2D $ shape2D
  --     initSourceSparse =
  --       computeInitialDistributionPowerMethodSparse phiFreqs rhoFreqs points
  --     harmonicsArraySparse =
  --       computeHarmonicsArraySparse
  --         numPointsRecon
  --         deltaRecon
  --         numPointsRecon
  --         deltaRecon
  --         phiFreqs
  --         rhoFreqs
  --         thetaFreqs
  --         scaleFreqs
  --         -- (2 * pi)
  --         (log (periodR2 * sqrt 2))
  --         64
  --         s
  -- ys <-
  --   computeContourSparse
  --     plan
  --     folderPath
  --     coefficients
  --     harmonicsArraySparse
  --     -- []
  --     (fromListUnboxed (Z :. L.length phiFreqs) phiFreqs)
  --     (fromListUnboxed (Z :. L.length rhoFreqs) rhoFreqs)
  --     64
  --     xs
  --     numIteration
  --     ""
  --     initSourceSparse
  -- let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
  --     idxVec = VU.fromList [(freqX, freqY) | freqX <- r2Freqs, freqY <- r2Freqs]
  --     eigenSourceSparse' =
  --       DFTArray numR2Freq numR2Freq phiFreqs rhoFreqs .
  --       L.map
  --         (L.foldl'
  --            (\vec ((x, y), z) ->
  --               VS.zipWith (+) vec .
  --               VU.convert .
  --               VU.map
  --                 (\(freqX, freqY) ->
  --                    z * (cis (-(freqX * x + freqY * y) * 2 * pi / periodR2))) $
  --               idxVec)
  --            (VS.replicate (numR2Freq ^ 2) 0) .
  --          L.zip xs) .
  --       L.transpose . L.map R.toList $
  --       ys
  --     eigenSourceSparse =
  --       convoluvePinhweelBasis coefficients harmonicsArray eigenSourceSparse'
  --     eigenSourceSparseMat =
  --       A.transpose . A.use . fromDFTArray . getDFTArrayVector $
  --       eigenSourceSparse
  -- eigenSourceSparseR2' <-
  --   computeFourierSeriesR2StreamAcc
  --     ptxs
  --     numR2Freq
  --     numPointsRecon
  --     ((2 * scaleFreq + 1) * (2 * thetaFreq + 1))
  --     periodR2
  --     deltaRecon
  --     numBatchR2
  --     eigenSourceSparseMat
  -- let eigenSourceSparseR2 =
  --       R.reshape
  --         (Z :. (2 * scaleFreq + 1) :. (2 * thetaFreq + 1) :. numPointsRecon :.
  --          numPointsRecon) $
  --       eigenSourceSparseR2'
  --     eigenSinkSparseR2 =
  --       timeReversalRepa
  --         (L.map fromIntegral . getListFromNumber $ (2 * thetaFreq + 1))
  --         eigenSourceSparseR2
  -- completionSparseR2 <-
  --   completionFieldRepa plan eigenSourceSparseR2 eigenSinkSparseR2
  -- plotImageRepa (folderPath </> "SourceSparsePinwheel.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
  --   VU.map sqrt .
  --   toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
  --   eigenSourceSparseR2
  -- plotImageRepa (folderPath </> "SinkSparsePinwheel.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
  --   VU.map sqrt .
  --   toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
  --   eigenSinkSparseR2
  -- plotImageRepa (folderPath </> "CompletionSparsePinwheel.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) .
  --   VU.map sqrt .
  --   toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
  --   completionSparseR2
  -- printCurrentTime "Discrete and sparse done"
  if scaleFreq == 0
    then let initialSourceDistribution =
               computeInitialDistributionPowerMethodPinwheelBasis'
                 numR2Freq
                 stdR2
                 periodR2
                 phiFreq
                 rhoFreq
                 points
         in computeContourPinwheelBasis'
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
    else let initialSourceDistribution =
               computeInitialDistributionPowerMethodPinwheelBasis
                 numR2Freq
                 stdR2
                 periodR2
                 phiFreq
                 rhoFreq
                 points
         in computeContourPinwheelBasis
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
