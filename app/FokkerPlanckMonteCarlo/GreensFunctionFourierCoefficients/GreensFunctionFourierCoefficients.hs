{-# LANGUAGE BangPatterns #-}
module GreensFunctionFourierCoefficients where

import           Data.Binary
import           Data.List               as L
import           FokkerPlanck.Histogram
import           FokkerPlanck.MonteCarlo
import           System.Directory
import           System.Environment
import           Utils.Time

main = do
  args@(thetaSigmaStr:scaleSigmaStr:maxScaleStr:deltaLogStr:taoStr:numTrailStr:maxTrailStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:histFilePath:numThreadStr:_) <-
    getArgs
  let thetaSigma = read thetaSigmaStr :: Double
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
      numThread = read numThreadStr :: Int
      folderPath = "output/app/GreensFunctionFourierCoefficients"
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      deltaLog = read deltaLogStr :: Double
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  initHist <-
    if flag
      then do printCurrentTime $
                "read Fourier coefficients data from" L.++ histFilePath
              decodeFile histFilePath
      else return
             (emptyHistogram
                [ L.length phiFreqs
                , L.length rhoFreqs
                , L.length thetaFreqs
                , L.length scaleFreqs
                ]
                0)
  runMonteCarloFourierCoefficients
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
    initHist
