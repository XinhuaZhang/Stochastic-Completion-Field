module STCR2Z2T0S0EndModal where

import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           Image.Transform         (normalizeValueRange)
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:cutoffRadiusEndPointStr:cutoffRadiusStr:reversalFactorStr:cStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      cutoffRadiusEndPoint = read cutoffRadiusEndPointStr :: Int
      cutoffRadius = read cutoffRadiusStr :: Int
      reversalFactor = read reversalFactorStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0EndModal.noindex"
      a = 20
      b = 7
      c = read cStr :: Int
      endPointFilePath =
        folderPath </>
        (printf
           "EndPoint_%d_%d_%d_%d_%d_%d_%.2f_%.2f_%d_%d_%d_%f.dat"
           numPoint
           (round thetaFreq :: Int)
           (round scaleFreq :: Int)
           (round maxScale :: Int)
           (round tao :: Int)
           cutoffRadiusEndPoint
           thetaSigma
           scaleSigma
           a
           b
           c
           reversalFactor)
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
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
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (cutoff cutoffRadius radialArr)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  arrR2Z2T0S0EndPoint <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (cutoff cutoffRadiusEndPoint radialArr)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  let a' = round $ (fromIntegral a) * (sqrt 2) / 2
      b' = round $ (fromIntegral b) * (sqrt 2) / 2
      c' = round $ (fromIntegral c) * (sqrt 2) / 2
      xs =
        ([R2S1RPPoint (i, i, 0, 0) | i <- [a',a' + b' .. c']] L.++
         [R2S1RPPoint (i, -i, 0, 0) | i <- [-a',-(a' + b') .. -c']] L.++
         [R2S1RPPoint (i, i, 0, 0) | i <- [-a',-(a' + b') .. -c']] L.++
         [R2S1RPPoint (i, -i, 0, 0) | i <- [a',a' + b' .. c']] L.++
         [R2S1RPPoint (i, 0, 0, 0) | i <- [a,a + b .. c]] L.++
         [R2S1RPPoint (0, i, 0, 0) | i <- [-a,-(a + b) .. -c]] L.++
         [R2S1RPPoint (i, 0, 0, 0) | i <- [-a,-(a + b) .. -c]] L.++
         [R2S1RPPoint (0, i, 0, 0) | i <- [a,a + b .. c]])
      bias = computeBiasR2T0S0' numPoint numPoint theta0Freqs scale0Freqs xs
      eigenVec =
        computeInitialEigenVectorR2T0S0
          numPoint
          numPoint
          theta0Freqs
          scale0Freqs
          thetaFreqs
          scaleFreqs
          xs
  powerMethodR2Z2T0S0EndModal
    plan
    folderPath
    numPoint
    numPoint
    numOrientation
    thetaFreqs
    theta0Freqs
    numScale
    scaleFreqs
    scale0Freqs
    arrR2Z2T0S0EndPoint
    arrR2Z2T0S0
    numIteration
    writeSourceFlag
    (printf
       "_%d_%d_%d_%d_%d_%d_%.2f_%.2f_%f"
       numPoint
       (round thetaFreq :: Int)
       (round scaleFreq :: Int)
       (round maxScale :: Int)
       (round tao :: Int)
       cutoffRadiusEndPoint
       thetaSigma
       scaleSigma
       reversalFactor)
    reversalFactor
    bias
    eigenVec
