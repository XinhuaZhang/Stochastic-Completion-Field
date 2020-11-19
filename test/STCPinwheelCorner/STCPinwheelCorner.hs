module STCPinwheelCorner where

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
import           FokkerPlanck
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
  args@(deviceIDsStr:numPointsStr:deltaStr:thresholdStr:numPointsReconStr:deltaReconStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:tauStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdR2Str:stdStr:numBatchR2Str:numBatchR2FreqsStr:numBatchOriStr:batchSizeStr:sStr:radiusStr:deltaTStr:weightStr:posStr:numThreadStr:_) <-
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
      folderPath = "output/test/STCPinwheelCorner"
      stdR2 = read stdR2Str :: Double
      std = read stdStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      numBatchOri = read numBatchOriStr :: Int
      batchSize = read batchSizeStr :: Int
      s = read sStr :: Double
      radius = read radiusStr :: Double
      deltaT = read deltaTStr :: Double
      weight = read weightStr :: Double
      pos = read posStr :: (Int, Int)
      periodEnv = periodR2 ^ 2 / 4
  -- removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  -- let deltaTheta = 2 * pi / fromIntegral numOrientation
  --     cornerDist =
  --       cornerDistribution
  --         delta
  --         initScale
  --         thetaSigma
  --         tau
  --         threshold
  --         [fromIntegral i * deltaTheta | i <- [0 .. numOrientation - 1]]
  --         weight
  --         pos
  --     idx = [(i, j) | i <- [-10 .. 10], j <- [-10 .. 10]]
  --     xs =
  --       L.map
  --         (\(i, j) ->
  --            let x = i * delta
  --                y = j * delta
  --             in ( (x, y)
  --                , L.sum $
  --                  cornerDistribution
  --                    delta
  --                    initScale
  --                    thetaSigma
  --                    tau
  --                    threshold
  --                    [ fromIntegral o * deltaTheta
  --                    | o <- [0 .. numOrientation - 1]
  --                    ]
  --                    weight
  --                    (x, y)))
  --         idx
  -- print . L.sortOn snd $ xs
  -- printf "%f\n" (L.sum cornerDist)
  -- print .
  --   L.sortOn snd .
  --   L.zip
  --     [fromIntegral i * deltaTheta / pi * 180 | i <- [0 .. numOrientation - 1]] $
  --   cornerDist
  -- plotImageRepa (folderPath </> "test.png") .
  --   ImageRepa 8 .
  --   fromListUnboxed (Z :. (1 :: Int) :. (21 :: Int) :. (21 :: Int)) .
  --   snd . L.unzip $
  --   xs
  -- let maxIdx =
  --       fst .
  --       L.maximumBy (\(_, a) (_, b) -> compare a b) .
  --       L.zip
  --         [ fromIntegral i * deltaTheta / pi * 180
  --         | i <- [0 .. numOrientation - 1]
  --         ] $
  --       cornerDist
  -- maxDist <- cornerDistribution' delta initScale thetaSigma tau threshold [fromIntegral i * deltaTheta | i <- [0 .. numOrientation - 1]] weight pos (maxIdx / 180 * pi)
  -- print maxIdx
  -- print maxDist
  flag <- doesFileExist histFilePath 
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  hist <-
    if flag
      then do
        printCurrentTime "read from files..."
        decodeFile histFilePath
      else sampleCartesianCorner
             histFilePath
             folderPath
             ptxs
             numPoints
             periodEnv
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
             stdR2
             weight
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
  let initDistSource =
        computeInitialDistributionFourierPinwheel
          numR2Freq
          periodR2
          periodEnv
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          initSource
      initDistSink =
        computeInitialDistributionFourierPinwheel
          numR2Freq
          periodR2
          periodEnv
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          initSink
  source <- convolve harmonicsArray initDistSource
  sink <- convolve harmonicsArray initDistSink
  sourceR2 <- plotFPArray plan (folderPath </> "Source.png") source
  sinkR2 <- plotFPArray plan (folderPath </> "Sink.png") sink
  plotRThetaDist
    (folderPath </> "Source_Theta.png")
    (folderPath </> "Source_R.png")
    numPointsRecon
    360
    90
    periodEnv
    pos
    sourceR2
  plotRThetaDist
    (folderPath </> "Sink_Theta.png")
    (folderPath </> "Sink_R.png")
    numPointsRecon
    360
    90
    periodEnv
    pos
    sinkR2
  sinkR2TR <-
    computeUnboxedP .
    timeReversalRepa [fromIntegral (-thetaFreq) .. fromIntegral thetaFreq] $
    sinkR2
  completionR2 <- completionFieldRepa plan sourceR2 sinkR2TR
  (sumP . sumS . R.map (\x -> magnitude x ** 2) . rotate4D2 $ completionR2) >>=
    plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) -- . R.map (\x -> log (x + 1))
