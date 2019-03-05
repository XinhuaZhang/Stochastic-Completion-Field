module STCR2Z2T0S0 where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary
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
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:initDistStr:histFilePath:numThreadStr:_) <-
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
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      initDist = read initDistStr :: [R2S1RPPoint]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numThread = read numThreadStr :: Int
      sourceDist = L.take 1 initDist
      sinkDist = L.drop 1 initDist
      folderPath = "output/test/STCR2Z2T0S0"
  flag <- doesFileExist histFilePath
  arrR2Z2T0S0 <-
    if flag
      then getNormalizedHistogramArr <$> decodeFile histFilePath
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
             histFilePath
             (emptyHistogram
                [ numPoint
                , numPoint
                , L.length scale0Freqs
                , L.length theta0Freqs
                , L.length scaleFreqs
                , L.length thetaFreqs
                ]
                0)
  createDirectoryIfMissing True folderPath
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  sourceDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      theta0Freqs
      scale0Freqs
      sourceDist
  sinkDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      theta0Freqs
      scale0Freqs
      sinkDist
  arrR2Z2T0S0F <- dftR2Z2T0S0 plan . makeFilterR2Z2T0S0 $ arrR2Z2T0S0
  -- Source field
  sourceArr <- convolveR2T0S0 plan arrR2Z2T0S0F sourceDistArr
  sourceR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D .
    r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs maxScale $
    sourceR2Z2
  plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- Sink field
  sinkArr <- convolveR2T0S0 plan arrR2Z2T0S0F sinkDistArr
  sinkR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sinkArr
  sinkField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D .
    r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs maxScale $
    sinkR2Z2
  plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  -- Completion Filed
  completionFiled <-
    timeReversalConvolveR2Z2 plan thetaFreqs scaleFreqs sourceR2Z2 sinkR2Z2
  completionFiledR2 <-
    R.sumP .
    R.sumS .
    R.map magnitude .
    rotate4D .
    rotate4D .
    r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs maxScale $
    completionFiled
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
