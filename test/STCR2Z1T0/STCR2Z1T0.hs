module STCR2Z1T0 where

import           Control.Monad             as M
import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
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
      folderPath = "output/test/STCR2Z1T0"
  flag <- doesFileExist histFilePath
  -- arrR2Z1T0 <-
  --   if pinwheelFlag
  --     then computeR2Z1T0Array numPoint numPoint alpha thetaFreqs theta0Freqs
  --     else if flag
  --            then getNormalizedHistogramArr <$> decodeFile histFilePath
  --            else solveMonteCarloR2Z1T0
  --                   numThread
  --                   numTrail
  --                   maxTrail
  --                   numPoint
  --                   numPoint
  --                   sigma
  --                   tao
  --                   1
  --                   len
  --                   theta0Freqs
  --                   thetaFreqs
  --                   histFilePath
  --                   (emptyHistogram
  --                      [ numPoint
  --                      , numPoint
  --                      , L.length theta0Freqs
  --                      , L.length thetaFreqs
  --                      ]
  --                      0)
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        solveMonteCarloR2Z1T0Radial
          numThread
          numTrail
          maxTrail
          numPoint
          numPoint
          sigma
          tao
          0
          theta0Freqs
          thetaFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length theta0Freqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z1T0 <-
    computeUnboxedP $
    computeR2Z1T0ArrayRadial
      radialArr
      numPoint
      numPoint
      1
      thetaFreqs
      theta0Freqs
  let arr3d =
        rotate3D . R.slice arrR2Z1T0 $
        (Z :. All :. (L.length theta0Freqs - 1) :. All :. All)
  createDirectoryIfMissing True folderPath
  MP.mapM_
    (\i ->
       plotImageRepaComplex
         (folderPath </> "GreensR2Z1T0_" L.++ show (i + 1) L.++ ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr3d $
       (Z :. All :. All :. i))
    [0 .. (L.length thetaFreqs) - 1]
  plan <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  sourceDistArr <-
    computeInitialDistributionR2T0 plan numPoint numPoint theta0Freqs sourceDist
  sinkDistArr <-
    computeInitialDistributionR2T0 plan numPoint numPoint theta0Freqs sinkDist
  arrR2Z1T0F <- (computeP . makeFilter2D $ arrR2Z1T0) >>= dftR2Z1T0 plan
  arrR2Z1T0TRF <-
    (computeP . makeFilter2D . timeReverseR2Z1T0 thetaFreqs theta0Freqs $
     arrR2Z1T0) >>=
    dftR2Z1T0 plan
  -- Source field
  sourceArr <- convolveR2T0 plan arrR2Z1T0F sourceDistArr
  -- let sourceR2Z1 = R.sumS . rotateR2Z1T0Array $ sourceArr
  --     sourceField =
  --       computeS .
  --       R.extend (Z :. (1 :: Int) :. All :. All) .
  --       R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
  --       sourceR2Z1
  sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sourceR2Z1
  plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- Sink field
  sinkArr' <- convolveR2T0 plan arrR2Z1T0F sinkDistArr
  let   sinkArr  = computeSinkFromSourceR2Z1T0 thetaFreqs theta0Freqs sinkArr'
  -- let sinkR2Z1 = R.sumS . rotateR2Z1T0Array $ sinkArr
  --     sinkField =
  --       computeS .
  --       R.extend (Z :. (1 :: Int) :. All :. All) .
  --       R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
  --       sinkR2Z1
  sinkR2Z1 <- R.sumP . rotateR2Z1T0Array $ sinkArr
  sinkField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sinkR2Z1
  plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  -- Completion Filed
  completionFiled <-
    convolveR2Z1 plan thetaFreqs sourceR2Z1 sinkR2Z1
  let completionFiledR2 =
        R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
        completionFiled
  -- completionFiledR2 <-
  --   R.sumP . R.map magnitude . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
  --   completionFiled
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
