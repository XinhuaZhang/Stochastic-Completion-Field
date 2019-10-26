module STCR2Z2T0S0Edge where

import           Control.Arrow
import           Control.Monad
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
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types


main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:thetaFreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:edgeFilePath:numNoisePointStr:scaleFactorStr:useFFTWWisdomFlagStr:fftwWisdomFileName:numThreadStr:_) <-
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
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Double
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      numNoisePoint = read numNoisePointStr :: Int
      scaleFactor = read scaleFactorStr :: Double
      useFFTWWisdomFlag = read useFFTWWisdomFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0Edge"
      fftwWisdomFilePath = folderPath </> fftwWisdomFileName
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
          thetaFreqs
          thetaFreqs
          scaleFreqs
          scaleFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scaleFreqs
             , L.length thetaFreqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (PinwheelHollow0 4)
      radialArr
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      thetaFreqs
      scaleFreqs
  plan <-
    makeR2Z2T0S0Plan emptyPlan useFFTWWisdomFlag fftwWisdomFilePath arrR2Z2T0S0
  xs <- parseEdgeFile edgeFilePath
  randomPonintSet <- generateRandomPointSet numNoisePoint numPoint numPoint
  let (centerX, centerY) =
        join (***) (\x -> div x . L.length $ xs) .
        L.foldl' (\(a, b) (R2S1RPPoint (c, d, _, _)) -> (a + c, b + d)) (0, 0) $
        xs
      ys =
        L.map
          (\(R2S1RPPoint (a, b, c, d)) ->
             (R2S1RPPoint
                ( round $ (fromIntegral $ a - centerX) / scaleFactor
                , round $ (fromIntegral $ b - centerX) / scaleFactor
                , c
                , d)))
          xs
      zs =
        L.map
          (\(R2S1RPPoint (a, b, c, d)) ->
             (R2S1RPPoint (a - center numPoint, b - center numPoint, c, d)))
          randomPonintSet
      points = ys L.++ zs
  let bias = computeBiasR2T0S0 numPoint numPoint thetaFreqs scaleFreqs points
      eigenVec =
        computeInitialEigenVectorR2T0S0
          numPoint
          numPoint
          thetaFreqs
          scaleFreqs
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
    numScale
    scaleFreqs
    maxScale
    arrR2Z2T0S0
    numIteration
    writeSourceFlag
    (printf
       "_%d_%d_%.2f_%.2f"
       (round maxScale :: Int)
       (round tao :: Int)
       thetaSigma
       scaleSigma)
    bias
    eigenVec
