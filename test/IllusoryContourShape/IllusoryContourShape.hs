{-# LANGUAGE BangPatterns #-}
module IllusoryContourShape where

import           Control.Monad               as M
import           Data.Array.Repa             as R
import           Data.Binary
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Generic         as VG
import           FokkerPlanck
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.Time

main = do
  args@(numPointStr:deltaXStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:deltaLogStr:taoStr:numTrailStr:maxTrailStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:shape2DStr:histFilePath:writeFlagStr:numIterationStr:numThreadStr:_) <-
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
      shape2D@(Points _ minDist _) = read shape2DStr :: Points Shape2D
      numThread = read numThreadStr :: Int
      folderPath =
        "output/test/IllusoryContourShape" </> takeBaseName histFilePath
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      deltaLog = read deltaLogStr :: Double
      writeFlag = read writeFlagStr :: Bool
      numIteration = read numIterationStr :: Int
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  hist <-
    if flag
      then do
        printCurrentTime $
          "read Fourier coefficients data from " L.++ histFilePath
        decodeFile histFilePath
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
  let !points =
        L.map (\(!x, !y) -> Point x y 0 1) . getShape2DIndexList . makeShape2D $
        shape2D
      !initSource =
        computeInitialDistributionPowerMethod
          numPoint
          numPoint
          phiFreqs
          rhoFreqs
          halfLogPeriod
          points
      !coefficients =
        getNormalizedHistogramArr hist :: R.Array U DIM4 (Complex Double)
      !thetaRHarmonics =
        computeThetaRHarmonics
          numOrientation
          numScale
          thetaFreqs
          scaleFreqs
          halfLogPeriod
      !bias = computeBias numPoint numPoint points
  plan <-
    makePlan
      emptyPlan
      numPoint
      numPoint
      (L.length thetaFreqs)
      (L.length scaleFreqs)
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
  completion <-
    computeContour
      plan
      folderPath
      writeFlag
      coefficients
      harmonicsArray
      bias
      numIteration
      initSource
  plotDFTArrayThetaR
    (folderPath </> (printf "Completion.png"))
    numPoint
    numPoint
    thetaRHarmonics
    completion
