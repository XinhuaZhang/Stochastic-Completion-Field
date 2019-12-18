module STCR2Z2T0S0 where

import           Control.Monad             as M
import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:initDistStr:histFilePath:alphaStr:radialFlagStr:numThreadStr:_) <-
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
      initDist = read initDistStr :: [R2S1RPPoint]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      radialFlag = read radialFlagStr :: Bool
      numThread = read numThreadStr :: Int
      sourceDist = L.take 1 initDist
      sinkDist = L.drop 1 initDist
      folderPath = "output/test/STCR2Z2T0S0"
      dirName = takeDirectory histFilePath
      fileName = takeFileName histFilePath
      newHistFilePath =
        if radialFlag
          then dirName </> ("Radial_" L.++ fileName)
          else histFilePath
      (R2S1RPPoint (_, _, _, s)) = L.head sourceDist
  createDirectoryIfMissing True folderPath
  arrR2Z2T0S0 <-
    if radialFlag
      then do
        flag <- doesFileExist newHistFilePath
        radialArr <-
          if flag
            then R.map magnitude . getNormalizedHistogramArr <$>
                 decodeFile newHistFilePath
            else do
              putStrLn
                "Couldn't find a Green's function data. Start simulation..."
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
                newHistFilePath
                (emptyHistogram
                   [ if maxScale >= 2
                        then round maxScale + 1
                        else (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
                   , L.length scale0Freqs
                   , L.length theta0Freqs
                   , L.length scaleFreqs
                   , L.length thetaFreqs
                   ]
                   0)
        computeUnboxedP $
          computeR2Z2T0S0ArrayRadial
            Pinwheel
            radialArr
            numPoint
            numPoint
            1
            maxScale
            thetaFreqs
            scaleFreqs
            theta0Freqs
            scale0Freqs
      else do
        flag <- doesFileExist newHistFilePath
        hist <-
          if flag
            then decodeFile newHistFilePath
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
        solveMonteCarloR2Z2T0S0
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
          newHistFilePath
          hist
  plan <- makeR2Z2T0S0Plan emptyPlan False "" arrR2Z2T0S0
  sourceDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      theta0Freqs
      scale0Freqs
      maxScale
      -- initDist
      sourceDist
  sinkDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      theta0Freqs
      scale0Freqs
      maxScale
      sinkDist
  arrR2Z2T0S0F <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ arrR2Z2T0S0
  -- Source field
  sourceArr <- convolveR2T0S0 plan arrR2Z2T0S0F sourceDistArr
  sourceR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D2 .
    r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs maxScale $
    sourceR2Z2
  plotImageRepaComplex (folderPath </> printf "Source%.0f.png" s) . ImageRepa 8 $
    sourceField
  -- Sink field
  sinkArr' <- convolveR2T0S0 plan arrR2Z2T0S0F sinkDistArr
  let sinkArr = computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sinkArr'
  sinkR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sinkArr
  sinkField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D2 .
    r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs maxScale $
    sinkR2Z2
  plotImageRepaComplex (folderPath </> printf "Sink%.0f.png" s) . ImageRepa 8 $
    sinkField
  -- Completion Filed
  completionFiled <- convolveR2Z2 plan sourceR2Z2 sinkR2Z2
  completionFiledR2 <-
    R.sumP .
    R.sumS .
    rotate4D2 .
    r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs maxScale $
    completionFiled
  plotImageRepaComplex (folderPath </> printf "Completion%.0f.png" s) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionMagnitude%.0f.png" s) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
