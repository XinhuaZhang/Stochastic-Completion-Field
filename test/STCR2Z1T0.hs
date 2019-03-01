module STCR2Z1T0 where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           STC.CompletionField
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:initDistStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      initDist = read initDistStr :: [R2S1RPPoint]
      numThread = read numThreadStr :: Int
      sourceDist = L.take 1 initDist
      sinkDist = L.drop 1 initDist
      folderPath = "output/test/STCR2Z1T0"
  arrR2Z1T0 <-
    solveMonteCarloR2Z1T0
      numThread
      numTrail
      maxTrail
      numPoint
      numPoint
      sigma
      tao
      len
      theta0Freqs
      thetaFreqs
      init
  createDirectoryIfMissing True folderPath
  plan <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  sourceDistArr <-
    computeInitialDistributionR2T0 plan numPoint numPoint theta0Freqs sourceDist
  sinkDistArr <-
    computeInitialDistributionR2T0 plan numPoint numPoint theta0Freqs sinkDist
  arrR2Z1T0F <- dftR2Z1T0 plan . makeFilterR2Z1T0 $ arrR2Z1T0
  -- Source field
  sourceArr <- convolveR2T0 plan arrR2Z1T0F sourceDistArr
  sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sourceR2Z1
  plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- Sink field
  sinkArr <- convolveR2T0 plan arrR2Z1T0F sinkDistArr
  sinkR2Z1 <- R.sumP . rotateR2Z1T0Array $ sinkArr
  sinkField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sinkR2Z1
  plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  -- Completion Filed
  completionFiled <-
    timeReversalConvolveR2Z1 plan thetaFreqs sourceR2Z1 sinkR2Z1
  completionFiledR2 <-
    R.sumP . R.map magnitude . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    completionFiled
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
