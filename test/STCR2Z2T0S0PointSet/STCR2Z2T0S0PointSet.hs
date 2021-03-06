module STCR2Z2T0S0PointSet where

import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC.PowerMethod
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import Text.Printf


main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:rStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
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
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0PointSet"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        solveMonteCarloR2Z2T0S0Radial
          numThread
          numTrail
          maxTrail
          numPoint
          numPoint
          thetaSigma
          scaleSigma
          maxScale
          tao
          theta0Freqs
          thetaFreqs
          scale0Freqs
          scaleFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (cutoff 8 radialArr)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  let n = 30
      m = round $ (fromIntegral n) * (sqrt 2) / 2
      a = 20
      b = 7
      c = 50
      a' = round $ (fromIntegral a) * (sqrt 2) / 2
      b' = round $ (fromIntegral b) * (sqrt 2) / 2
      c' = round $ (fromIntegral c) * (sqrt 2) / 2
      r = read rStr :: Double
      numTheta = 8
      deltaTheta = (2 * pi) / numTheta
      xs = -- [R2S1RPPoint (i,0,0,1) | i <- [-30,-20..30]]
        -- ((L.map
        --     (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
        --     [ (r * cos (k * deltaTheta) + 0, r * sin (k * deltaTheta) + 0)
        --     | k <- [0 .. numTheta-1]
        --     ]) -- L.++
         -- (L.map
         --    (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
         --    [ ((r + 2) * cos (k * deltaTheta), (r + 2) * sin (k * deltaTheta))
         --    | k <- [0 .. numTheta - 1]
         --    ]) L.++
         -- (L.map
         --    (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
         --    [ ( (r + 5) * cos (k * deltaTheta)
         --      , (r + 5) * sin (k * deltaTheta))
         --    | k <- [0 .. numTheta - 1]
         --    ])
               -- L.++ [R2S1RPPoint (3, 5, 0, 1)]
         -- )
               ([R2S1RPPoint (i, i, 0, 0) | i <- [a',a' + b' .. c']] L.++
                [R2S1RPPoint (i, -i, 0, 0) | i <- [-a',-(a' + b') .. -c']] L.++
                [R2S1RPPoint (i, i, 0, 0) | i <- [-a',-(a' + b') .. -c']] L.++
                [R2S1RPPoint (i, -i, 0, 0) | i <- [a',a' + b' .. c']] L.++
                [R2S1RPPoint (i, 0, 0, 0) | i <- [a,a + b .. c]] L.++
                [R2S1RPPoint (0, i, 0, 0) | i <- [-a,-(a + b) .. -c]] L.++
                [R2S1RPPoint (i, 0, 0, 0) | i <- [-a,-(a + b) .. -c]] L.++
                [R2S1RPPoint (0, i, 0, 0) | i <- [a,a + b .. c]] -- L.++
                -- [R2S1RPPoint (3, 5, 0, 1)]
               )
                -- [ r2s1rppoint (n, 0, 0, 1)
                -- , R2S1RPPoint (0, n, 0, 1)
                -- , R2S1RPPoint (-n, 0, 0, 1)
                -- , R2S1RPPoint (0, -n, 0, 1)
                -- , R2S1RPPoint (m, m, 0, 1)
                -- , R2S1RPPoint (-m, m, 0, 1)
                -- , R2S1RPPoint (m, -m, 0, 1)
                -- , R2S1RPPoint (-m, -m, 0, 1)
                -- ]
  let bias = computeBiasR2T0S0 numPoint numPoint theta0Freqs scale0Freqs xs
      eigenVec =
        computeInitialEigenVectorR2T0S0
          numPoint
          numPoint
          theta0Freqs
          scale0Freqs
          thetaFreqs
          scaleFreqs
          xs
  powerMethodR2Z2T0S0Reversal
    plan
    folderPath
    numPoint
    numPoint
    numOrientation
    thetaFreqs
    theta0Freqs
    numScale
    scaleFreqs
    scale0Freqs
    arrR2Z2T0S0
    numIteration
    writeSourceFlag
    -- (printf "_%d" (round r :: Int))
    (printf
       "_%d_%d_%d_%d_%d_%.2f_%.2f"
       numPoint
       (round thetaFreq :: Int)
       (round scaleFreq :: Int)
       (round maxScale :: Int)
       (round tao :: Int)
       thetaSigma
       scaleSigma)
    0.5
    0.001
    bias
    eigenVec
