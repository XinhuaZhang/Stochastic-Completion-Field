module STCR2Z2T0S0Image where

import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types


main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:imagePath:cutoffRStr:reversalFactorStr:numThreadStr:_) <-
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
      cutoffR = read cutoffRStr :: Int
      reversalFactor = read reversalFactorStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0Image"
  imgRepa@(ImageRepa _ img') <- readImageRepa imagePath False
  let (Z :. _ :. cols :. rows) = extent img'
      img =
        computeUnboxedS . R.traverse img' id $ \f idx@(Z :. k :. i :. j) ->
          if (sqrt . fromIntegral $ (i - div cols 2) ^ 2 + (j - div rows 2) ^ 2) >
             48
            then 0
            else f idx
  createDirectoryIfMissing True folderPath
  plotImageRepa (folderPath </> "input.png") . ImageRepa 8 $ img
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
          cols
          rows
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
             [ (round . sqrt . fromIntegral $ 2 * (div (max rows cols) 2) ^ 2)
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (cutoff cutoffR radialArr)
      cols
      rows
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  powerMethodR2Z2T0S0Reversal
    plan
    folderPath
    cols
    rows
    numOrientation
    thetaFreqs
    theta0Freqs
    numScale
    scaleFreqs
    scale0Freqs
    arrR2Z2T0S0
    numIteration
    writeSourceFlag
    (printf
       "_%d_%d_%.2f_%.2f"
       (round maxScale :: Int)
       (round tao :: Int)
       thetaSigma
       scaleSigma)
    0.5
    reversalFactor
    (R.traverse
       img
       (const
          (Z :. (L.length theta0Freqs) :. (L.length scale0Freqs) :. cols :. rows)) $ \f (Z :. _ :. _ :. i :. j) ->
       if f (Z :. (0 :: Int) :. i :. j) > 0
         then 1
         else 0)
    (R.traverse
       img
       (const
          (Z :. (L.length thetaFreqs) :. (L.length scaleFreqs) :.
           (L.length theta0Freqs) :.
           (L.length scale0Freqs) :.
           cols :.
           rows)) $ \f (Z :. k :. l :. _ :. _ :. i :. j) ->
       if k == div (L.length thetaFreqs) 2 && l == div (L.length scaleFreqs) 2
         then 1 / (fromIntegral $ (L.length theta0Freqs * L.length scale0Freqs)) :+
              0
         else 0)
