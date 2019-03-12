module STCR2Z2T0S0Image where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           Filter.Pinwheel
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC.CompletionField
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:alphaStr:pinwheelFlagStr:imagePath:numThreadStr:_) <-
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
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numThread = read numThreadStr :: Int
      pinwheelParams =
        PinwheelParams numPoint numPoint alpha maxScale theta0Freqs scale0Freqs
      folderPath = "output/test/STCR2Z2T0S0Image"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  arrR2Z2T0S0 <-
    if pinwheelFlag
      then error
             "Using pinwheels to construct the Green's function has not been implemented yet."
      else if flag
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
  
  (ImageRepa _ img) <- readImageRepa imagePath False
  plan0 <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  (plan1, imgF) <- makeImagePlan plan0 img
  (plan, filterF, filterPIF) <- pinwheelFilter plan1 pinwheelParams
  arrR2Z2T0S0F <- dftR2Z2T0S0 plan . makeFilterR2Z2T0S0 $ arrR2Z2T0S0
  let initialDistF =
        computeS $
        frequencyDomainMultiply
          filterF
          (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  -- Source field
  sourceArr <- convolveR2T0S0 plan arrR2Z2T0S0F initialDistF
  sourceR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    sourceR2Z2
  plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- Completion Filed
  completionFiled <-
    timeReversalConvolveR2Z2 plan thetaFreqs scaleFreqs sourceR2Z2 sourceR2Z2 
  completionFiledR2 <-
    R.sumP .
    R.sumS .
    R.map magnitude .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    completionFiled
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
