module IllusoryContourEdgeCorner where
import           Control.Arrow
import           Control.Monad                  as M
import           Control.Monad.Parallel         as MP
-- import qualified Data.Array.Accelerate          as A
-- import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           FokkerPlanck
-- import           Foreign.CUDA.Driver            as CUDA
-- import           FourierMethod.BlockMatrixAcc
-- import           FourierMethod.FourierSeries2D
import           Image.Edge
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
import Utils.Distribution
import FourierPinwheel.GaussianEnvelopePinwheel
import Filter.Utils

{-# INLINE dropPoints #-}
dropPoints :: [a] -> [a]
dropPoints [] = []
dropPoints (x:[]) = []
dropPoints (x:y:[]) = [x]
dropPoints (x:y:xs) = x : dropPoints xs

{-# INLINE findMinMax #-}
findMinMax :: (Ord a) => [(a, a)] -> (a -> a -> a) -> a
findMinMax [] _ = error "findMinMax: empty list"
findMinMax ((x, y):[]) f = f x y
findMinMax ((x, y):xs) f = f (f x y) (findMinMax xs f)

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:tauStr:numTrailsStr:deltaTStr:poissonWeightStr:numR2FreqStr:periodR2Str:deltaFreqStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initScaleStr:histFilePath:stdR2Str:stdThetaStr:stdRStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:writeFlagStr:numIterationStr:radiusStr:edgeFilePath1:edgeFilePath2:edgeFilePath3:scaleFactorStr:numNoisePointStr:noisePointRangeStr:numThreadStr:_) <-
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
      poissonWeight = read poissonWeightStr :: Double
      numR2Freq = read numR2FreqStr :: Int
      periodR2 = read periodR2Str :: Double
      deltaFreq = read deltaFreqStr :: Double
      phiFreq = read phiFreqsStr :: Int
      phiFreqs = getListFromNumber phiFreq
      rhoFreq = read rhoFreqsStr :: Int
      rhoFreqs = L.map (* deltaFreq) . getListFromNumber' $ rhoFreq
      thetaFreq = read thetaFreqsStr :: Int
      thetaFreqs = getListFromNumber thetaFreq
      scaleFreq = read scaleFreqsStr :: Int
      scaleFreqs = L.map (* deltaFreq) . getListFromNumber' $ scaleFreq
      initScale = read initScaleStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/IllusoryContourEdgeCorner"
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
      radius = read radiusStr :: Double
      scaleFactor = read scaleFactorStr :: Double
      numNoisePoint = read numNoisePointStr :: Int
      noisePointRange = read noisePointRangeStr :: (Double, Double)
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  hist <-
    if flag
      then do
        printCurrentTime $ "read coefficients from file \n" L.++ histFilePath
        decodeFile histFilePath
      else do
        printCurrentTime "Start computing coefficients..."
        sampleCartesianCorner
          histFilePath
          folderPath
          numThread
          numPoints
          delta
          numOrientation
          initScale
          thetaSigma
          tau
          threshold
          s
          phiFreqs
          rhoFreqs
          thetaFreqs
          scaleFreqs
          stdR2
          poissonWeight
      -- else do
      --   printCurrentTime "Start computing coefficients..."
      --   initialise []
      --   devs <- M.mapM device deviceIDs
      --   ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
      --   ptxs <- M.mapM createTargetFromContext ctxs
      --   sampleCartesian
      --     histFilePath
      --     folderPath
      --     ptxs
      --     numPoints
      --     delta
      --     numOrientation
      --     initScale
      --     thetaSigma
      --     tau
      --     threshold
      --     s
      --     phiFreqs
      --     rhoFreqs
      --     thetaFreqs
      --     scaleFreqs
      --     stdR2
  printCurrentTime "Done"
  printCurrentTime "Start Convloution.."
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPointsRecon
      numPointsRecon
      numR2Freq
      thetaFreq
      scaleFreq
      phiFreq
      rhoFreq
  printCurrentTime "Compute DFT Plan done."
  let coefficients = getNormalizedHistogramArr hist
      asteriskGaussianVec =
        gaussianPinwheel
          numR2Freq
          periodR2
          stdR2
          10
          thetaFreqs
          scaleFreqs
          stdTheta
          stdR
  harmonicsArray <-
    createHarmonics
      numR2Freq
      phiFreq
      rhoFreq
      thetaFreq
      scaleFreq
      deltaFreq
      (-s)
      periodR2
      coefficients
  M.mapM_
    (\filePath -> do
       xs' <- parseEdgeFile' filePath
       -- let xs = dropPoints $ dropPoints $ dropPoints $ dropPoints $ dropPoints $ L.zip x'  y'
       let newFolderPath = folderPath </> takeBaseName filePath
       createDirectoryIfMissing True newFolderPath
       let xs = dropPoints $ dropPoints $ dropPoints $ dropPoints $ xs'
           (centerX, centerY) =
             join (***) (\x -> x / (fromIntegral . L.length $ xs)) .
             L.foldl' (\(a, b) (c, d) -> (a + c, b + d)) (0, 0) $
             xs
           ys =
             L.map
               (\(a, b) ->
                  ((a - centerX) * scaleFactor, (centerY - b) * scaleFactor))
               xs
       print (centerX, centerY)
       let points = L.map (\(a, b) -> Point b a 0 1) ys
       (bias1, dftBias1) <-
         computeBiasFourierPinwheelFull
           plan
           numR2Freq
           thetaFreq
           scaleFreq
           (-s)
           periodR2
           radius
           stdTheta
           stdR
           stdR2
           asteriskGaussianVec
           points
       let initDist1 =
             computeInitialDistributionPowerMethodFourierPinwheelFull
               numR2Freq
               phiFreq
               rhoFreq
               thetaFreq
               scaleFreq
               bias1
       plotFPArray'
         plan
         (newFolderPath </> "Source.png")
         initDist1
         (1 / scaleFactor)
         (-(centerX - 128), (centerY - 128))
       (bias, dftBias) <-
         computeBiasFourierPinwheelFull
           plan
           numR2Freq
           thetaFreq
           scaleFreq
           (-s)
           periodR2
           radius
           stdTheta
           stdR
           stdR2
           asteriskGaussianVec
           points
       let initDist =
             computeInitialDistributionPowerMethodFourierPinwheelFull
               numR2Freq
               phiFreq
               rhoFreq
               thetaFreq
               scaleFreq
               bias
       plotFPArray plan (newFolderPath </> "SourceNoise.png") initDist
       computeContourFourierPinwheelRecenterRescale
         plan
         newFolderPath
         writeFlag
         harmonicsArray
         dftBias
         numIteration
         numBatchR2
         numPointsRecon
         deltaRecon
         periodR2
         initDist
         ""
         (1 / scaleFactor)
         (-(centerX - 128), (centerY - 128)))
    [edgeFilePath2, edgeFilePath3]
