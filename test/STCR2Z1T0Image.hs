module STCR2Z1T0Image where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
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
import Data.Vector.Generic as VG
import Data.Vector.Storable as VS
import System.Random

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:histFilePath:alphaStr:pinwheelFlagStr:imagePath:numIterationStr:numThreadStr:_) <-
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
      numIteration = read numIterationStr :: Int
      numThread = read numThreadStr :: Int
      pinwheelParams =
        PinwheelParams numPoint numPoint alpha (exp 1) theta0Freqs [0]
      folderPath = "output/test/STCR2Z1T0Image"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  arrR2Z1T0 <-
    if pinwheelFlag
      then error
             "Using pinwheels to construct the Green's function has not been implemented yet."
           --computeR2Z1T0Array numPoint numPoint alpha thetaFreqs theta0Freqs
      else if flag
               -- getNormalizedHistogramArr <$> decodeFile histFilePath
             then do
               h <-  decodeFile histFilePath
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
                 histFilePath
                 h
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
  let arr3d =
        rotate3D . R.slice arrR2Z1T0 $
        (Z :. All :. (L.length theta0Freqs - 1) :. All :. All)
  M.mapM_
    (\i ->
       plotImageRepaComplex
         (folderPath </> "GreensR2Z1T0_" L.++ show (i + 1) L.++ ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr3d $
       (Z :. All :. All :. i))
    [0 .. (L.length thetaFreqs) - 1]
  (ImageRepa _ img) <- readImageRepa imagePath False
  -- let img =
  --       computeS $
  --       fromFunction
  --         (Z :. (1 :: Int) :. numPoint :. numPoint)
  --         (\(Z :. _ :. i :. j) ->
  --            if (i == (div numPoint 2) -- || i == (div numPoint 2) - 8 || i == (div numPoint 2) + 8
  --               ) &&
  --               ((abs (j - (div numPoint 2 - 12)) < 6) ||
  --                (abs (j - (div numPoint 2 + 12)) < 6)
  --               )
  --              then 255
  --              else 0)
  plotImageRepa (folderPath </> "input.png") (ImageRepa 8 img)
  plan0 <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  (plan1, imgF) <- makeImagePlan plan0 img
  (plan, filterF, filterPIF) <- pinwheelFilter plan1 pinwheelParams
  arrR2Z1T0F <- dftR2Z1T0 plan . makeFilterR2Z1T0 $ arrR2Z1T0
  convolvedImg <-
    convolvePinwheel plan filterF (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  convolvedImgPI <-
    convolvePinwheel
      plan
      filterPIF
      (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  let initialDist =
        computeUnboxedS . R.slice convolvedImg $
        (Z :. All :. (0 :: Int) :. All :. All)
      initialDistPI =
        computeUnboxedS . R.slice convolvedImgPI $
        (Z :. All :. (0 :: Int) :. All :. All)
  xs1Source <-
    M.replicateM ((L.length theta0Freqs) * numPoint * numPoint) randomIO
  xs2Source <-
    M.replicateM ((L.length theta0Freqs) * numPoint * numPoint) randomIO
  let eigenVecSource =
        fromListUnboxed (Z :. (L.length theta0Freqs) :. numPoint :. numPoint) $
        L.zipWith (:+) xs1Source xs2Source
  arrR2Z1Source <-
    M.foldM
      (\input n ->
         eigenVectorR2Z1Source
           plan
           folderPath
           numOrientation
           thetaFreqs
           arrR2Z1T0F
           n
           (computeS . R.slice convolvedImg $
            (Z :. All :. (0 :: Int) :. All :. All))
           input)
      eigenVecSource
      [1 .. numIteration]
  xs1Sink <-
    M.replicateM ((L.length theta0Freqs) * numPoint * numPoint) randomIO
  xs2Sink <-
    M.replicateM ((L.length theta0Freqs) * numPoint * numPoint) randomIO
  let eigenVecSink =
        fromListUnboxed (Z :. (L.length theta0Freqs) :. numPoint :. numPoint) $
        L.zipWith (:+) xs1Sink xs2Sink
  arrR2Z1Sink <-
    M.foldM
      (\input n ->
         eigenVectorR2Z1Sink
           plan
           folderPath
           numOrientation
           thetaFreqs
           arrR2Z1T0F
           n
           (computeS . R.slice convolvedImgPI $
            (Z :. All :. (0 :: Int) :. All :. All))
           input)
      eigenVecSink
      [1 .. numIteration]
  completionFieldR2Z1
    plan
    folderPath
    numOrientation
    thetaFreqs
    arrR2Z1Source
    arrR2Z1Sink
  -- let sourceInitialDistF =
  --       computeS $
  --       R.slice
  --         (frequencyDomainMultiply
  --            filterPIF
  --            (R.slice imgF (Z :. (0 :: Int) :. All :. All)))
  --         (Z :. All :. (0 :: Int) :. All :. All)
      -- sinkInitialDistF =
      --   computeS $
      --   R.slice
      --     (frequencyDomainMultiply
      --        filterPIF
      --        (R.slice imgF (Z :. (0 :: Int) :. All :. All)))
      --     (Z :. All :. (0 :: Int) :. All :. All)
  -- convolvedImg <-
  --   convolvePinwheel plan filterF (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  -- convolvedImgPI <-
  --   convolvePinwheel
  --     plan
  --     filterPIF
  --     (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  -- let initialDist =
  --       R.slice convolvedImg $
  --       (Z :. All :. (0 :: Int) :. All :. All)
  --     initialDistPI =
  --       R.slice convolvedImgPI $
  --       (Z :. All :. (0 :: Int) :. All :. All)
  -- sourceInitialDistF <-
  --   fmap
  --     (fromUnboxed (Z :. (L.length theta0Freqs) :. numPoint :. numPoint) .
  --      VG.convert) .
  --   dftExecute
  --     plan
  --     (DFTPlanID DFT1DG [(L.length theta0Freqs), numPoint, numPoint] [1, 2]) .
  --   VG.convert . toUnboxed . computeS $ -- initialDistPI
  --   R.zipWith (+) initialDist initialDistPI
  -- sinkInitialDistF <-
  --   fmap
  --     (fromUnboxed (Z :. (L.length theta0Freqs) :. numPoint :. numPoint) .
  --      VG.convert) .
  --   dftExecute
  --     plan
  --     (DFTPlanID DFT1DG [(L.length theta0Freqs), numPoint, numPoint] [1, 2]) .
  --   VG.convert . toUnboxed . computeS $
  --   initialDistPI
  -- -- Source field
  -- sourceArr <- convolveR2T0 plan arrR2Z1T0F sourceInitialDistF
  -- sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceArr
  -- sourceField <-
  --   fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
  --   R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
  --   sourceR2Z1
  -- plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- -- -- Sink field
  -- -- sinkArr <- convolveR2T0 plan arrR2Z1T0F sinkInitialDistF
  -- -- sinkR2Z1 <- R.sumP . rotateR2Z1T0Array $ sinkArr
  -- -- sinkField <-
  -- --   fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
  -- --   R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
  -- --   sinkR2Z1
  -- -- plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  -- -- Completion Filed
  -- completionFiled <-
  --   timeReversalConvolveR2Z1 plan thetaFreqs sourceR2Z1 sourceR2Z1
  -- completionFiledR2 <-
  --   R.sumP . R.map magnitude . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
  --   completionFiled
  -- plotImageRepa (folderPath </> "Completion.png") .
  --   ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
  --   completionFiledR2
