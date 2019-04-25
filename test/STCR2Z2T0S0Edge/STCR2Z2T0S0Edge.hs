module STCR2Z2T0S0Edge where

import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.Edge
import           Image.IO
import           STC.CompletionField
import           STC.PointSet
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types


main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:edgeFilePath:scaleFactorStr:numThreadStr:_) <-
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
      scaleFactor = read scaleFactorStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0Edge"
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
      radialArr
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  xs <- parseEdgeFile edgeFilePath
  randomPonintSet <- generateRandomPointSet 40 numPoint numPoint
  let ys =
        L.map
          (\(R2S1RPPoint (a, b, c, d)) ->
             (R2S1RPPoint
                ( round $ (fromIntegral a - 150) / scaleFactor
                , round $ (fromIntegral b - 150) / scaleFactor
                , c
                , d)))
          xs
      zs =
        L.map
          (\(R2S1RPPoint (a, b, c, d)) ->
             (R2S1RPPoint (a - center numPoint, b - center numPoint, c, d)))
          randomPonintSet
      points = ys L.++ zs
  let bias = computeBiasR2T0S0 numPoint numPoint theta0Freqs scale0Freqs points
      eigenVec =
        computeInitialEigenVectorR2T0S0
          numPoint
          numPoint
          theta0Freqs
          scale0Freqs
          thetaFreqs
          scaleFreqs
          points
  powerMethodR2Z2T0S0
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
    arrR2Z2T0S0
    numIteration
    writeSourceFlag
    -- (printf "_%d" (round scaleFactor :: Int))
    (printf "_%d_%d_%.2f_%.2f" (round maxScale :: Int) (round tao :: Int) thetaSigma scaleSigma)
    0.5
    bias
    eigenVec
