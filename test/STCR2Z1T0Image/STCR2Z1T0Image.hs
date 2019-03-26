module STCR2Z1T0Image where

import           Control.Monad                as M
import           Control.Monad.Parallel       as MP
import           Data.Array.Repa              as R
import           Data.Binary                  (decodeFile)
import           Data.Complex
import           Data.List                    as L
import           Data.Vector.Generic          as VG
import           Data.Vector.Storable         as VS
import           DFT.Plan
import           Filter.Gaussian
import           Filter.Pinwheel
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           Image.Transform
import           STC.CompletionField
import           STC.OrientationScaleAnalysis
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initialScaleStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:histFilePath:alphaStr:pinwheelFlagStr:imagePath:numIterationStr:writeSourceFlagStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      initialScale = read initialScaleStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      numThread = read numThreadStr :: Int
      pinwheelParams =
        PinwheelParams numPoint numPoint alpha (exp 1) theta0Freqs [0]
      folderPath = "output/test/STCR2Z1T0Image"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  arrR2Z1T0' <-
    if pinwheelFlag
      then error
             "Using pinwheels to construct the Green's function has not been implemented yet."
          --computeR2Z1T0Array numPoint numPoint alpha thetaFreqs theta0Freqs
      else if flag
              -- getNormalizedHistogramArr <$> decodeFile histFilePath
              {-h <-  decodeFile histFilePath
              solveMonteCarloR2Z1T0
                numThread
                numTrail
                maxTrail
                numPoint
                numPoint
                sigma
                tao
                initialScale
                len
                theta0Freqs
                thetaFreqs
                histFilePath
                h-}
             then do
               getNormalizedHistogramArr <$> decodeFile histFilePath
             else do print "Couldn't find a Green's function data. Start simulation...\n"
                     solveMonteCarloR2Z1T0
                       numThread
                       numTrail
                       maxTrail
                       numPoint
                       numPoint
                       sigma
                       tao
                       initialScale
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
  let arrR2Z1T0
        -- computeUnboxedS .
        -- pad [numPoint, numPoint, L.length theta0Freqs, L.length thetaFreqs] 0 .
        -- downsample [2, 2, 1, 1] $
       = arrR2Z1T0'
  -- -- let arr3d =
  -- --       rotate3D . R.slice arrR2Z1T0 $
  -- --       (Z :. All :. (L.length theta0Freqs - 1) :. All :. All)
  -- -- MP.mapM_
  -- --    (\i ->
  -- --       plotImageRepaComplex
  -- --         (folderPath </> "GreensR2Z1T0_" L.++ show (i + 1) L.++ ".png") .
  -- --       ImageRepa 8 .
  -- --       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr3d $
  -- --       (Z :. All :. All :. i))
  -- --    [0 .. (L.length thetaFreqs) - 1]
  (ImageRepa _ img') <- readImageRepa imagePath False
  plan0 <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  (plan1, imgF) <- makeImagePlan plan0 img'
  (plan2, filterF, filterPIF) <- pinwheelFilter plan1 pinwheelParams
  (plan, gaussianFilterF) <-
    gaussian2DFilter plan2 (Gaussian2DParams 1 numPoint numPoint)
  -- -- img <- convolveGaussian2D plan gaussianFilterF . R.map (:+ 0) $ img'
  -- -- imgVec <- dropPixel (2 / 3) . toUnboxed . computeS . filterImage $ img'
  let img = R.map (:+ 0) img' -- = . fromUnboxed (extent img') $ imgVec
  plotImageRepa
    (folderPath </> "input.png") -- (ImageRepa 8 img)
    (ImageRepa 8 . computeS . R.map magnitude $ img)
  convolvedImg' <-
    convolvePinwheel plan filterPIF (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  let convolvedImg =
        computeS . R.slice convolvedImg' $
        (Z :. All :. (0 :: Int) :. All :. All)
  let (initialDistSource, initialDistSink) =
        analyzeOrientation
          numOrientation
          theta0Freqs
          convolvedImg
  powerMethod2
    plan
    folderPath
    numPoint
    numOrientation
    thetaFreqs
    theta0Freqs
    arrR2Z1T0
    numIteration
    writeSourceFlag
    initialDistSource
