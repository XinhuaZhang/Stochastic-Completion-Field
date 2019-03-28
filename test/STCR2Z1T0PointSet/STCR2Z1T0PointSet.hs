module STCR2Z1T0PointSet where

import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC.CompletionField
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types


main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initialScaleStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:histFilePath:pinwheelFlagStr:numIterationStr:writeSourceFlagStr:numThreadStr:_) <-
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
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z1T0PointSet"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
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
      1
      thetaFreqs
      theta0Freqs
  plan <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  let n = 30
      m = round $ (fromIntegral n) * (sqrt 2) / 2
      a = 30
      b = -30
      c = 30
      a' = round $ (fromIntegral a) * (sqrt 2) / 2
      b' = round $ (fromIntegral b) * (sqrt 2) / 2
      c' = round $ (fromIntegral c) * (sqrt 2) / 2
  initialDistSource <-
    computeInitialDistributionR2T0
      plan
      numPoint
      numPoint
      theta0Freqs
    -- [R2S1RPPoint (i, 0, 0, 1) | i <- [-20,-18 .. (-10)]]
      ([R2S1RPPoint (i, i, 180, 1) | i <- [a',a' + b' .. c']] L.++
       [R2S1RPPoint (i, -i, 0, 1) | i <- [-a',-(a' + b') .. -c']] L.++
       [R2S1RPPoint (i, i, 180, 1) | i <- [-a',-(a' + b') .. -c']] L.++
       [R2S1RPPoint (i, -i, 0, 1) | i <- [a',a' + b' .. c']] L.++
       [R2S1RPPoint (i, 0, 180, 1) | i <- [a,a + b .. c]] L.++
       [R2S1RPPoint (0, i, 0, 1) | i <- [-a,-(a + b) .. -c]] L.++
       [R2S1RPPoint (i, 0, 180, 1) | i <- [-a,-(a + b) .. -c]] L.++
       [R2S1RPPoint (0, i, 0, 1) | i <- [a,a + b .. c]] L.++ [R2S1RPPoint (3, 5, 0, 0)])
       -- [ R2S1RPPoint (n, 0, 0, 1)
       -- , R2S1RPPoint (0, n, 0, 1)
       -- , R2S1RPPoint (-n, 0, 0, 1)
       -- , R2S1RPPoint (0, -n, 0, 1)
       -- , R2S1RPPoint (m, m, 0, 1)
       -- , R2S1RPPoint (-m, m, 0, 1)
       -- , R2S1RPPoint (m, -m, 0, 1)
       -- , R2S1RPPoint (-m, -m, 0, 1)
       -- ]
  powerMethod1
    plan
    folderPath
    numPoint
    numPoint
    numOrientation
    thetaFreqs
    theta0Freqs
    arrR2Z1T0
    numIteration
    writeSourceFlag
    ""
    initialDistSource
