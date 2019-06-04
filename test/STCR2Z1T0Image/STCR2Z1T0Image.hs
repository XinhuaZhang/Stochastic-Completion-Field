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
import           STC
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
      folderPath = "output/test/STCR2Z1T0Image"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  (ImageRepa _ img'') <- readImageRepa imagePath False
  let img'
        -- resize25D (numPoint, numPoint) (0, 255) .
        -- R.backpermute (Z :. nf :. numPoint :. numPoint) id .
        -- upsample [4, 4, 1] . downsample [4, 4, 1] $
       = img''
      (Z :. nf :. _ :. _) = extent img''
      (Z :. _ :. cols :. rows) = extent img'
      pinwheelParams = PinwheelParams rows cols alpha (exp 1) theta0Freqs [0]
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
          1
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
      1.5
      thetaFreqs
      theta0Freqs
  plan <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  -- (plan1, imgF) <- makeImagePlan plan0 . computeS $ img'
  -- (plan2, filterF, filterPIF) <- pinwheelFilter plan1 pinwheelParams
  -- (plan, gaussianFilterF) <-
  --   gaussian2DFilter plan2 (Gaussian2DParams 1 rows cols)
  -- -- img <- convolveGaussian2D plan gaussianFilterF . R.map (:+ 0) $ img'
  -- -- imgVec <- dropPixel (2 / 3) . toUnboxed . computeS . filterImage $ img'
  let img = R.map (:+ 0) img' -- = . fromUnboxed (extent img') $ imgVec
  plotImageRepa
    (folderPath </> "input.png") -- (ImageRepa 8 img)
    (ImageRepa 8 . computeS . R.map magnitude $ img)
  -- convolvedImg' <-
  --   convolvePinwheel
  --     plan
  --     filterPIF
  --     (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  -- let convolvedImg =
  --       computeS . R.slice convolvedImg' $
  --       (Z :. All :. (0 :: Int) :. All :. All)
  -- let (initialDistSource, initialDistSink) =
  --       analyzeOrientation numOrientation theta0Freqs convolvedImg
  powerMethod1
    plan
    folderPath
    cols
    rows
    numOrientation
    thetaFreqs
    theta0Freqs
    arrR2Z1T0
    numIteration
    writeSourceFlag
    ""
    0.1
    (R.traverse img (const (Z :. (L.length theta0Freqs) :. cols :. rows)) $ \f (Z :. _ :. i :. j) ->
       if magnitude (f (Z :. (0 :: Int) :. i :. j)) > 0
         then 1
         else 0)
    (R.traverse
       img
       (const
          (Z :. (L.length thetaFreqs) :. (L.length theta0Freqs) :. cols :. rows)) $ \f (Z :. k :. _ :. i :. j) ->
       if k == div (L.length thetaFreqs) 2
         then f (Z :. (0 :: Int) :. i :. j) /
              (fromIntegral $ L.length theta0Freqs)
         else 0)
