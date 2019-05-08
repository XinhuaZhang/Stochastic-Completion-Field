module STCR2Z2T0S0ReversalCornerPointSet where

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
import           Text.Printf
import           Types


main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoDecayStr:taoReversalStr:taoCornerStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:rStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      taoDecay = read taoDecayStr :: Double
      taoReversal = read taoReversalStr :: Double
      taoCorner = read taoCornerStr :: Double
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
      folderPath = "output/test/STCR2Z2T0S0ReversalCornerPointSet"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        solveMonteCarloR2Z2T0S0ReversalCornerRadial
          numThread
          numTrail
          maxTrail
          numPoint
          numPoint
          thetaSigma
          scaleSigma
          maxScale
          taoDecay
          taoReversal
          taoCorner
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
      radialArr
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
      a = 10
      b = -10
      c = 10
      a' = round $ (fromIntegral a) * (sqrt 2) / 2
      b' = round $ (fromIntegral b) * (sqrt 2) / 2
      c' = round $ (fromIntegral c) * (sqrt 2) / 2
      r = read rStr :: Double
      rr = round r :: Int
      shift = 15
      numTheta = 45
      deltaTheta = 1 * pi / numTheta
      xs =
        -- [ R2S1RPPoint (rr, shift, 0, 1)
        -- , R2S1RPPoint (rr, -shift, 0, 1)
        -- , R2S1RPPoint (-rr, shift, 0, 1)
        -- , R2S1RPPoint (-rr, -shift, 0, 1)
        -- , R2S1RPPoint (shift, rr, 0, 1)
        -- , R2S1RPPoint (-shift, rr, 0, 1)
        -- , R2S1RPPoint (shift, -rr, 0, 1)
        -- , R2S1RPPoint (-shift, -rr, 0, 1)
        -- ]
        -- [R2S1RPPoint (i, -i, 0, 1) | i <- [-2,0 .. 20]] L.++
        [R2S1RPPoint (-i, -i, 0, 1) | i <- [-10,-8 .. 10]]
        -- ((L.map
        --     (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
        --     [ (r * cos (k * deltaTheta) + 0, r * sin (k * deltaTheta) + 0)
        --     | k <- [0 .. numTheta]
        --     ]) -- L.++
        --  -- (L.map
        --  --    (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
        --  --    [ ((r + 2) * cos (k * deltaTheta), (r + 2) * sin (k * deltaTheta))
        --  --    | k <- [0 .. numTheta - 1]
        --  --    ]) L.++
        --  -- (L.map
        --  --    (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
        --  --    [ ( (r + 5) * cos (k * deltaTheta)
        --  --      , (r + 5) * sin (k * deltaTheta))
        --  --    | k <- [0 .. numTheta - 1]
        --  --    ])
        --        -- L.++ [R2S1RPPoint (3, 5, 0, 1)]
        --  )
               -- ([R2S1RPPoint (i, i, 0, 1) | i <- [a',a' + b' .. c']] L.++
               --  [R2S1RPPoint (i, -i, 0, 1) | i <- [-a',-(a' + b') .. -c']] L.++
               --  [R2S1RPPoint (i, i, 0, 1) | i <- [-a',-(a' + b') .. -c']] L.++
               --  [R2S1RPPoint (i, -i, 0, 1) | i <- [a',a' + b' .. c']] L.++
               --  [R2S1RPPoint (i, 0, 0, 1) | i <- [a,a + b .. c]] L.++
               --  [R2S1RPPoint (0, i, 0, 1) | i <- [-a,-(a + b) .. -c]] L.++
               --  [R2S1RPPoint (i, 0, 0, 1) | i <- [-a,-(a + b) .. -c]] L.++
               --  [R2S1RPPoint (0, i, 0, 1) | i <- [a,a + b .. c]] L.++
               --  [R2S1RPPoint (3, 5, 0, 1)])
                -- [ R2S1RPPoint (n, 0, 0, 1)
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
  powerMethodR2Z2T0S0
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
       "_%d_%d_%d_%d_%.2f_%.2f"
       (round maxScale :: Int)
       (round taoDecay :: Int)
       (round taoReversal :: Int)
       (round taoCorner :: Int)
       thetaSigma
       scaleSigma)
    0.5
    bias
    eigenVec
