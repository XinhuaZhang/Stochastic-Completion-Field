{-# LANGUAGE FlexibleContexts #-}
module R2Z2T0S0ToR2S1RPT0S0 where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.Time

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histPath:alphaStr:pinwheelFlagStr:numThreadStr:_) <-
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
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/R2Z2T0S0ToR2S1RPT0S0"
  flag <- doesFileExist histPath
  -- arrR2Z2T0S0 <-
  --   -- if flag
  --   --   then getNormalizedHistogramArr <$> decodeFile histPath
  --   --   else
  --     solveMonteCarloR2Z2T0S0
  --            numThread
  --            numTrail
  --            maxTrail
  --            numPoint
  --            numPoint
  --            thetaSigma
  --            scaleSigma
  --            maxScale
  --            tao
  --            len
  --            theta0Freqs
  --            thetaFreqs
  --            scale0Freqs
  --            scaleFreqs
  --            histPath
  --            (emptyHistogram
  --               [ numPoint
  --               , numPoint
  --               , L.length scale0Freqs
  --               , L.length theta0Freqs
  --               , L.length scaleFreqs
  --               , L.length thetaFreqs
  --               ]
  --               0)
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$> decodeFile histPath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        printCurrentTime
        solveMonteCarloR2Z2T0S0Radial
          numThread
          numTrail
          maxTrail
          numPoint
          numPoint
          thetaSigma
          scaleSigma
          maxScale
          tao
          theta0Freqs
          thetaFreqs
          scale0Freqs
          scaleFreqs
          ""
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  printCurrentTime
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      radialArr
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  arrR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ arrR2Z2T0S0
      -- arr =
      --   r2z2t0s0Tor2s1rpt0s0
      --     numOrientation
      --     thetaFreqs
      --     numScale
      --     scaleFreqs
      --     arrR2Z2T0S0
      -- arr4d =
      --   R.slice arr $
      --   (Z :. All :. All :. (L.length theta0Freqs - 1) :.
      --    (L.length scale0Freqs - 1) :.
      --    All :.
      --    All)
  let arr4d = r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs arrR2Z2
  createDirectoryIfMissing True folderPath
  printCurrentTime
  MP.mapM_
    (\(i, j) ->
       plotImageRepaComplex
         (folderPath </> show (i + 1) L.++ "_" L.++ show (j + 1) L.++ ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr4d $
       (Z :. i :. j :. All :. All))
    [(i, j) | i <- [0 .. numOrientation - 1], j <- [0 .. numScale - 1]]
  arr2D <- R.sumP . R.sumS . rotate4D . rotate4D $ arr4d
  plotImageRepaComplex (folderPath </> "sum.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    arr2D
