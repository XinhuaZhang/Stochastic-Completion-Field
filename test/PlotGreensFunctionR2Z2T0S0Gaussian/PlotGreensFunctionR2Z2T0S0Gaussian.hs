{-# LANGUAGE FlexibleContexts #-}
module PlotGreensFunctionR2Z2T0S0Gaussian where

import           Control.Monad.Parallel  as MP
import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.List               as L
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFileName:alphaStr:gaussianSigmaStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      thetaFreq = read thetaFreqsStr :: Double
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      gaussianSigma = read gaussianSigmaStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/PlotGreensFunctionR2Z2T0S0Gaussian"
      histPath = folderPath </> histFileName
  createDirectoryIfMissing True (folderPath </> "GreensR2Z2T0S0Gaussian")
  flag <- doesFileExist histPath
  hist <-
    if flag
      then decodeFile histPath
      else return $
           emptyHistogram
             [ numPoint
             , numPoint
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0
  arr <-
    solveMonteCarloR2Z2T0S0Gaussian
      numThread
      numTrail
      maxTrail
      numPoint
      numPoint
      thetaSigma
      scaleSigma
      maxScale
      tao
      len
      theta0Freqs
      thetaFreqs
      scale0Freqs
      scaleFreqs
      gaussianSigma
      histPath
      hist
  let arr4d =
        R.slice
          arr
          (Z :. All :. All :. (L.length theta0Freqs - 1) :.
           (L.length scale0Freqs - 1) :.
           All :.
           All)
  MP.mapM_
    (\(i, j) ->
       plotImageRepaComplex
         (folderPath </> "GreensR2Z2T0S0Gaussian" </> show (i + 1) L.++ "_" L.++
          show (j + 1) L.++
          ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr4d $
       (Z :. i :. j :. All :. All))
    [ (i, j)
    | i <- [0 .. (L.length thetaFreqs) - 1]
    , j <- [0 .. (L.length scaleFreqs) - 1]
    ]
