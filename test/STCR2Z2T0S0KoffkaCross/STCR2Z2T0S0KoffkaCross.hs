module STCR2Z2T0S0KoffkaCross where

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
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoDecayStr:taoReversalStr:taoCornerStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:dStr:wStr:sigmaStr:thetaStr:numThreadStr:_) <-
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
      d = read dStr :: Int
      w = read wStr :: Int
      sigma = read sigmaStr :: Double
      theta = read thetaStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0KoffkaCross"
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
  let xs =
        [ R2S1RPPoint (d, w, theta, 1)
        , R2S1RPPoint (d, -w, theta, 1)
        , R2S1RPPoint (w, d, theta + 90, 1)
        , R2S1RPPoint (-w, d, theta + 90, 1)
        , R2S1RPPoint (-d, w, theta + 180, 1)
        , R2S1RPPoint (-d, -w, theta + 180, 1)
        , R2S1RPPoint (w, -d, theta + 270, 1)
        , R2S1RPPoint (-w, -d, theta + 270, 1)
        ]
      bias =
        -- computeBiasR2T0S0Gaussian
        --   numPoint
        --   numPoint
        --   theta0Freqs
        --   scale0Freqs
        --   90
        --   sigma
        --   xs
        computeS $
        R.zipWith
          (+)
          (computeBiasR2T0S0Gaussian
             numPoint
             numPoint
             theta0Freqs
             scale0Freqs
             90
             sigma
             xs)
          (computeBiasR2T0S0Gaussian
             numPoint
             numPoint
             theta0Freqs
             scale0Freqs
             (-90)
             sigma
             xs)
      eigenVec =
        computeInitialEigenVectorR2T0S0
          numPoint
          numPoint
          theta0Freqs
          scale0Freqs
          thetaFreqs
          scaleFreqs
          xs
  powerMethodR2Z2T0S0Bias
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
       "_%d_%d_%d_%d_%d_%.2f_%.2f_%d_%d"
       numPoint
       (round maxScale :: Int)
       (round taoDecay :: Int)
       (round taoReversal :: Int)
       (round taoCorner :: Int)
       thetaSigma
       scaleSigma
       d
       w)
    0.5
    bias
    eigenVec
