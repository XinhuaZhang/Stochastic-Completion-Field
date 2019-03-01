{-# LANGUAGE FlexibleContexts #-}
module R2Z2T0S0ToR2S1RPT0S0 where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.List                 as L
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histPath:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
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
      numThread = read numThreadStr :: Int
      folderPath = "output/test/R2Z2T0S0ToR2S1RPT0S0"
  flag <- doesFileExist histPath
  arrR2Z2T0S0 <-
    if flag
      then getNormalizedHistogramArr <$> decodeFile histPath
      else solveMonteCarloR2Z2T0S0
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
             histPath
             (emptyHistogram
                [ numPoint
                , numPoint
                , L.length scale0Freqs
                , L.length theta0Freqs
                , L.length scaleFreqs
                , L.length thetaFreqs
                ]
                0)
  let arr =
        r2z2t0s0Tor2s1rpt0s0
          numOrientation
          thetaFreqs
          numScale
          scaleFreqs
          maxScale
          arrR2Z2T0S0
      arr4d =
        R.slice arr $
        (Z :. All :. All :. (L.length theta0Freqs - 1) :.
         (L.length scale0Freqs - 1) :.
         All :.
         All)
  createDirectoryIfMissing True folderPath
  MP.mapM_
    (\(i, j) ->
       plotImageRepaComplex
         (folderPath </> show (i + 1) L.++ "_" L.++ show (j + 1) L.++ ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr4d $
       (Z :. i :. j :. All :. All))
    [(i, j) | i <- [0 .. numOrientation - 1], j <- [0 .. numScale - 1]]
