module OrientationAnalysis where

import           Control.Monad                as M
import           Data.Array.Repa              as R
import           Data.Binary                  (decodeFile)
import           Data.Complex
import           Data.List                    as L
import           DFT.Plan
import           Filter.Pinwheel
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC.CompletionField          (makeImagePlan, makeR2Z1T0Plan)
import           STC.CompletionField
import           STC.OrientationScaleAnalysis
import           System.Directory
import           System.Environment
import           System.FilePath

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:histFilePath:alphaStr:pinwheelFlagStr:imagePath:idxStr:numOrientationSampleStr:numThreadStr:_) <-
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
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      idx = read idxStr :: (Int, Int)
      numOrientationSample = read numOrientationSampleStr :: Int
      numThread = read numThreadStr :: Int
      pinwheelParams =
        PinwheelParams numPoint numPoint alpha (exp 1) theta0Freqs [0]
      folderPath = "output/test/OrientationAnalysis"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  arrR2Z1T0 <-
    if pinwheelFlag
      then error
             "Using pinwheels to construct the Green's function has not been implemented yet."
           --computeR2Z1T0Array numPoint numPoint alpha thetaFreqs theta0Freqs
      else if flag
             then do
               getNormalizedHistogramArr <$> decodeFile histFilePath
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
                    ""
                    (emptyHistogram
                       [ numPoint
                       , numPoint
                       , L.length theta0Freqs
                       , L.length thetaFreqs
                       ]
                       0)
  (ImageRepa _ img) <- readImageRepa imagePath False
  plotImageRepa (folderPath </> "input.png") (ImageRepa 8 img)
  plan0 <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  (plan1, imgF) <- makeImagePlan plan0 img
  (plan, filterF, filterPIF) <- pinwheelFilter plan1 pinwheelParams
  arrR2Z1T0F <- dftR2Z1T0 plan . makeFilterR2Z1T0 $ arrR2Z1T0
  convolvedImg <-
    fmap (\x -> R.slice x (Z :. All :. (0 :: Int) :. All :. All)) $
    convolvePinwheel plan filterF (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  plotMagnitudeOrientationSource
    plan
    folderPath
    numOrientationSample
    numOrientation
    thetaFreqs
    arrR2Z1T0F
    convolvedImg
    idx
  plotMagnitudeOrientation
    folderPath
    numOrientationSample
    numOrientation
    thetaFreqs
    convolvedImg
    idx
