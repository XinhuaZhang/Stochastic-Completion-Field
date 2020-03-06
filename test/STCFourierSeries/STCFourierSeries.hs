{-# LANGUAGE BangPatterns #-}
module STCFourierSeries where

import           Control.DeepSeq
import           Control.Monad                      as M
import qualified Data.Array.Accelerate              as A
import qualified Data.Array.Accelerate.Data.Complex as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.IArray
import           Data.Array.Repa                    as R
import           Data.Binary
import           Data.Complex
import           Data.List                          as L
import           Data.Maybe
import           Data.Vector.Generic                as VG
import           FokkerPlanck.BrownianMotion
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.Histogram
import           FokkerPlanck.MonteCarlo
import           Foreign.CUDA.Driver                as CUDA
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.Time

main = do
  args@(gpuIDStr:numPointStr:deltaXStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:deltaLogStr:taoStr:numTrailStr:maxTrailStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdStr:numThreadStr:_) <-
    getArgs
  let gpuID = read gpuIDStr :: [Int]
      numPoint = read numPointStr :: Int
      deltaX = read deltaXStr :: Double
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      phiFreq = read phiFreqsStr :: Double
      phiFreqs = [-phiFreq .. phiFreq]
      rhoFreq = read rhoFreqsStr :: Double
      rhoFreqs = [-rhoFreq .. rhoFreq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Double
      scaleFreqs = [-scaleFreq .. scaleFreq]
      initDist = read initDistStr :: [(Int, Int, Double, Double)]
      initScale = read initScaleStr :: Double
      initPoints = L.map (\(x, y, t, s) -> Point x y t s) initDist
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCFourierSeries"
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      deltaLog = read deltaLogStr :: Double
      std = read stdStr :: Double
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  hist <-
    if flag
      then do
        printCurrentTime $ "read data from " L.++ histFilePath
        decodeFile histFilePath
      else runMonteCarloFourierCoefficientsGPU
             gpuID
             numThread
             numTrail
             maxTrail
             thetaSigma
             scaleSigma
             maxScale
             tao
             phiFreqs
             rhoFreqs
             thetaFreqs
             scaleFreqs
             deltaLog
             initScale
             histFilePath
             (emptyHistogram
                [ L.length phiFreqs
                , L.length rhoFreqs
                , L.length thetaFreqs
                , L.length scaleFreqs
                ]
                0)
  let !initSource =
        computeInitialDistribution'
          numPoint
          numPoint
          phiFreqs
          rhoFreqs
          -- thetaFreqs
          halfLogPeriod
          [L.head initPoints]
      !initSink =
        computeInitialDistribution'
          numPoint
          numPoint
          phiFreqs
          rhoFreqs
          -- thetaFreqs
          halfLogPeriod
          [L.last initPoints]
      !coefficients =
        normalizeFreqArr' std phiFreqs rhoFreqs .
        getNormalizedHistogramArr $
        hist
       -- = getNormalizedHistogramArr $ hist :: R.Array U DIM4 (Complex Double)
      !thetaRHarmonics =
        computeThetaRHarmonics
          numOrientation
          numScale
          thetaFreqs
          scaleFreqs
          halfLogPeriod
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPoint
      numPoint
      (L.length thetaFreqs)
      (L.length scaleFreqs)
  printCurrentTime "harmonicsArray"
  harmonicsArray <-
    dftHarmonicsArray
      plan
      numPoint
      deltaX
      numPoint
      deltaX
      phiFreqs
      rhoFreqs
      thetaFreqs
      scaleFreqs
      halfLogPeriod
      (fromIntegral numPoint)
  -- let harmonicsArray' =
  --       computeHarmonicsArray
  --         numPoint
  --         deltaX
  --         numPoint
  --         deltaX
  --         phiFreqs
  --         rhoFreqs
  --         thetaFreqs
  --         scaleFreqs
  --         halfLogPeriod
  --         maxScale
  --         -- (fromIntegral numPoint)
  -- harmonicsArrayGPU <-
  --   dftHarmonicsArrayGPU
  --     plan
  --     numPoint
  --     deltaX
  --     numPoint
  --     deltaX
  --     phiFreqs
  --     rhoFreqs
  --     thetaFreqs
  --     scaleFreqs
  --     halfLogPeriod
  --     maxScale
  -- initialise []
  -- dev <- device . L.head $ gpuID
  -- ctx <- CUDA.create dev []
  -- ptx <- createTargetFromContext ctx
  -- registerPinnedAllocatorWith ptx
  --Source
  printCurrentTime "Source"
  sourceArr <- convolve' Source plan coefficients harmonicsArray initSource
  -- sourceArr' <- convolve'' Source plan coefficients harmonicsArray' initSource
  -- printCurrentTime "CPU"
  -- sourceArr <-
  --   convolveGPU
  --     ptx
  --     Source
  --     plan
  --     (A.use .
  --      A.fromList
  --        (A.Z A.:. (L.length scaleFreqs) A.:. (L.length thetaFreqs) A.:.
  --         (L.length rhoFreqs) A.:.
  --         (L.length phiFreqs)) .
  --      R.toList $
  --      coefficients)
  --     (A.use .
  --      A.fromList
  --        (A.Z A.:. (L.length scaleFreqs) A.:. (L.length thetaFreqs) A.:.
  --         (L.length rhoFreqs) A.:.
  --         (L.length phiFreqs)) .
  --      R.toList $
  --      coefficients)
  --     harmonicsArrayGPU
  --     initSource
  -- printCurrentTime "GPU"
  plotDFTArrayPower
    (folderPath </> (printf "SourcePower.png"))
    numPoint
    numPoint $
    sourceArr
  plotDFTArrayThetaR
    (folderPath </> "Source.png")
    numPoint
    numPoint
    thetaRHarmonics
    sourceArr
  -- plotDFTArrayThetaR
  --   (folderPath </> "Source_test.png")
  --   numPoint
  --   numPoint
  --   thetaRHarmonics
  --   sourceArr'
  --Sink
  printCurrentTime "Sink"
  sinkArr <- convolve' Sink plan coefficients harmonicsArray initSink
  plotDFTArrayPower (folderPath </> (printf "SinkPower.png")) numPoint numPoint $
    sinkArr
  plotDFTArrayThetaR
    (folderPath </> "Sink.png")
    numPoint
    numPoint
    thetaRHarmonics
    sinkArr
  printCurrentTime "Completion"
  completionArr <- completionField' plan sourceArr sinkArr
  plotDFTArrayThetaR
    (folderPath </> "Completion.png")
    numPoint
    numPoint
    thetaRHarmonics
    completionArr
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower.png"))
    numPoint
    numPoint
    completionArr
  let a =
        computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $
        sourceArr
      b =
        computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $ sinkArr
  plotImageRepa (folderPath </> "CompletionR2S1.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) .
    VG.convert . L.foldl1' (VG.zipWith (+)) $
    L.zipWith (VG.zipWith (\x y -> magnitude x * magnitude y)) a b
  printCurrentTime ""
  -- M.mapM_
  --   (\s -> do
  --      let !arr =
  --            getNormalizedHistogramArr hist :: R.Array U DIM4 (Complex Double)
  --          initArr =
  --            traverse4
  --              (fromListUnboxed (Z :. L.length phiFreqs) phiFreqs)
  --              (fromListUnboxed (Z :. L.length rhoFreqs) rhoFreqs)
  --              (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
  --              (fromListUnboxed (Z :. L.length scaleFreqs) scaleFreqs)
  --              (\_ _ _ _ -> extent arr) $ \fPhi fRho fTheta fScale (Z :. scale :. theta :. rho :. phi) ->
  --              cis $
  --              (-pi * log s) * fRho (Z :. rho) / halfLogPeriod -
  --              (fPhi (Z :. phi) * t0 * pi / 180)
  --      printCurrentTime ""
  --      source <- computeUnboxedP $ arr *^ initArr
  --      plotImageRepaComplex (folderPath </> (printf "Source_%.2f.png" s)) .
  --        ImageRepa 8 .
  --        fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) .
  --        L.foldl1' (VU.zipWith (+)) .
  --        computeFourierSeriesThetaR
  --          numOrientation
  --          numScale
  --          thetaFreqs
  --          scaleFreqs
  --          halfLogPeriod .
  --        computeFourierSeriesR2
  --          numPoint
  --          deltaX
  --          numPoint
  --          deltaX
  --          phiFreqs
  --          rhoFreqs
  --          thetaFreqs
  --          scaleFreqs
  --          halfLogPeriod
  --          maxScale
  --          harmonicsArray $
  --        source)
  --   [1,1.25 .. s0]
