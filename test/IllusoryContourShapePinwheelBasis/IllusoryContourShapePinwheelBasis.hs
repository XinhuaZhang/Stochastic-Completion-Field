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
import Utils.Distribution
import FourierPinwheel.GaussianEnvelopePinwheel
import Filter.Utils

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:tauStr:numTrailsStr:deltaTStr:poissonWeightStr:numR2FreqStr:periodR2Str:deltaFreqStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:histFilePathCorner:stdR2Str:stdThetaStr:stdRStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:writeFlagStr:numIterationStr:shape2DStr:radiusStr:numThreadStr:_) <-
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
      initDist = read initDistStr :: [(Double, Double, Double, Double)]
      initPoints = L.map (\(x, y, t, s) -> Point x y t s) initDist
      initSource = [L.head initPoints]
      initSink = [L.last initPoints]
      numThread = read numThreadStr :: Int
      folderPath = "output/test/IllusoryContourShapePinwheelBasis"
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
      shape2Ds = read shape2DStr :: [Points Shape2D]
      radius = read radiusStr :: Double
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  flagCorner <- doesFileExist histFilePathCorner
  hist <-
    case (getShape . L.head $ shape2Ds) of
      Circle _ _ ->
        if flag
          then do
            printCurrentTime $ "read coefficients from file \n" L.++ histFilePath
            decodeFile histFilePath
          else do
            printCurrentTime "Start computing coefficients..."
            initialise []
            devs <- M.mapM device deviceIDs
            ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
            ptxs <- M.mapM createTargetFromContext ctxs
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
              phiFreqs
              rhoFreqs
              thetaFreqs
              scaleFreqs
              stdR2
      KoffkaCross _ _ ->
        if flagCorner
          then do
            printCurrentTime "read coefficients from file"
            decodeFile histFilePathCorner
          else do
            printCurrentTime "Start computing coefficients..."
            initialise []
            devs <- M.mapM device deviceIDs
            ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
            ptxs <- M.mapM createTargetFromContext ctxs
            sampleCartesianCorner
              histFilePathCorner
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
              phiFreqs
              rhoFreqs
              thetaFreqs
              scaleFreqs
              stdR2
              poissonWeight
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
        case getShape (L.head shape2Ds) of
          Circle _ _ ->
            gaussianPinwheel
              numR2Freq
              periodR2
              stdR2
              10
              thetaFreqs
              scaleFreqs
              stdTheta
              stdR
          KoffkaCross _ _ ->
            gaussianPinwheelDiscontinuity
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
    (\shape2D -> do
       let points =
             L.map (\(x, y) -> Point (x) (y) 0 1) . -- L.map (translate (0.7, 0)) . L.map (rotateNdilate 30 1) .
             getShape2DIndexList' . makeShape2D $
             shape2D
       printCurrentTime (show points)
       (bias, dftBias) <-
         case getShape shape2D of
           Circle _ _ ->
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
           KoffkaCross _ _ ->
             computeBiasFourierPinwheelKoffkaCorss
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
       initialise []
       devs <- M.mapM device deviceIDs
       ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
       ptxs <- M.mapM createTargetFromContext ctxs
       arr <-plotFPArrayAcc
               ptxs
               (folderPath </> (show . getShape $ shape2D) L.++ "_Source.png" )
               numPointsRecon
               deltaRecon
               periodR2
               numBatchR2
               initDist
       plotRThetaDist (folderPath </> ("Theta_Source.png" )) (folderPath </> ("R_Source.png" ))  numPointsRecon 180 90 (exp (2*pi)) (36, 2) arr
       computeContourFourierPinwheel
         plan
         folderPath
         writeFlag
         harmonicsArray
         dftBias
         numIteration
         numBatchR2
         numPointsRecon
         deltaRecon
         periodR2
         initDist
         (show . getShape $ shape2D)
         deviceIDs)
    shape2Ds
