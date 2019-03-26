module STCR2Z1T0_1 where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           Image.Transform
import           STC.CompletionField
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:initDistStr:histFilePath:alphaStr:pinwheelFlagStr:numThreadStr:_) <-
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
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numThread = read numThreadStr :: Int
      sourceDist = L.take 1 initDist
      sinkDist = L.drop 1 initDist
      folderPath = "output/test/STCR2Z1T0_1"
  flag <- doesFileExist histFilePath
  arrR2Z1T0' <-
    if pinwheelFlag
      then computeR2Z1T0Array numPoint numPoint alpha thetaFreqs theta0Freqs
      else if flag
             then getNormalizedHistogramArr <$> decodeFile histFilePath
             else solveMonteCarloR2Z1T0
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
                    histFilePath
                    (emptyHistogram
                       [ numPoint
                       , numPoint
                       , L.length theta0Freqs
                       , L.length thetaFreqs
                       ]
                       0)
  let arrR2Z1T0 =
        computeUnboxedS .
        pad [numPoint, numPoint, L.length theta0Freqs, L.length thetaFreqs] 0 .
        downsample [2, 2, 1, 1] $
        arrR2Z1T0'
  createDirectoryIfMissing True folderPath
  plan <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  sourceDistArr <-
    computeInitialDistributionR2Z1T0
      plan
      numPoint
      numPoint
      thetaFreqs
      theta0Freqs
      sourceDist
  sinkDistArr <-
    computeInitialDistributionR2Z1T0
      plan
      numPoint
      numPoint
      thetaFreqs
      theta0Freqs
      sinkDist
  arrR2Z1T0F <- dftR2Z1T0 plan . makeFilterR2Z1T0 $ arrR2Z1T0
  arrR2Z1T0TRF <-
    dftR2Z1T0 plan . makeFilterR2Z1T0 . timeReverseR2Z1T0 thetaFreqs theta0Freqs $
    arrR2Z1T0
  -- Source field
  sourceArr <- convolveR2Z1T0 plan arrR2Z1T0F sourceDistArr
  sourceR2 <- R.sumP . R.sumS . rotate4D . rotate4D $ sourceArr
  let sourceField =
        computeS . R.extend (Z :. (1 :: Int) :. All :. All) $ sourceR2
  plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- Sink field
  sinkArr <- convolveR2Z1T0 plan arrR2Z1T0TRF sinkDistArr
  sinkR2 <- R.sumP . R.sumS . rotate4D . rotate4D $ sinkArr
  let sinkField = computeS . R.extend (Z :. (1 :: Int) :. All :. All) $ sinkR2
  plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  -- Completion Filed
  completionFiled <- convolveR2Z1' plan thetaFreqs theta0Freqs sourceArr sinkArr
  completionFiledR2 <- R.sumP . R.sumS . rotate4D . rotate4D $ completionFiled
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
