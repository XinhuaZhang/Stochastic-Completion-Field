{-# LANGUAGE BangPatterns #-}
module STCFourierSeries where

import           Control.Monad               as M
import           Data.Array.Repa             as R
import           Data.Binary
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Generic         as VG
import           FokkerPlanck.BrownianMotion
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.Histogram
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           STC.Convolution
import           STC.InitialDistribution
import           STC.Plan
import           STC.Point
import           STC.CompletionField
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.Time

main = do
  args@(numPointStr:deltaXStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:deltaLogStr:taoStr:numTrailStr:maxTrailStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:histFilePath:numThreadStr:_) <-
    getArgs
  let numPoint = read numPointStr :: Int
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
      initPoints = L.map (\(x, y, t, s) -> Point x y t s) initDist
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCFourierSeries"
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      deltaLog = read deltaLogStr :: Double
      -- (_, _, t0, s0) = initDist
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  hist <-
    if flag
      then decodeFile histFilePath
      else runMonteCarloFourierCoefficients
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
             histFilePath
             (emptyHistogram
                [ L.length phiFreqs
                , L.length rhoFreqs
                , L.length thetaFreqs
                , L.length scaleFreqs
                ]
                0)
  let !initSource =
        computeInitialDistribution
          numPoint
          numPoint
          phiFreqs
          rhoFreqs
          halfLogPeriod
          [L.head initPoints]
      !initSink =
        computeInitialDistribution
          numPoint
          numPoint
          phiFreqs
          rhoFreqs
          halfLogPeriod
          [L.last initPoints]
      !coefficients =
        getNormalizedHistogramArr hist :: R.Array U DIM4 (Complex Double)
      !thetaRHarmonics =
        computeThetaRHarmonics
          numOrientation
          numScale
          thetaFreqs
          scaleFreqs
          halfLogPeriod
  plan <-
    makePlan
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
      maxScale
  --Source
  printCurrentTime "Source"
  sourceArr <- convolve Source plan coefficients harmonicsArray initSource
  -- plotImageRepaComplex (folderPath </> (printf "Source.png")) .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) .
  --   VG.convert .
  --   L.foldl1' (VG.zipWith (+)) .
  --   computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $
  --   sourceArr
  --Sink
  printCurrentTime "Sink"
  sinkArr <- convolve Sink plan coefficients harmonicsArray initSink
  -- plotImageRepaComplex (folderPath </> (printf "Sink.png")) .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) .
  --   VG.convert .
  --   L.foldl1' (VG.zipWith (+)) .
  --   computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $
  --   sinkArr
  printCurrentTime "Completion"
  completionArr <- completionField plan sourceArr sinkArr
  let !convertedArr =
        L.foldl1' (VG.zipWith (+)) .
        computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $
        completionArr
  plotImageRepaComplex (folderPath </> (printf "Completion.png")) .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) . VG.convert $
    convertedArr
  plotImageRepa (folderPath </> (printf "CompletionMag.png")) .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) .
    VG.convert . VG.map magnitude $
    convertedArr
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
