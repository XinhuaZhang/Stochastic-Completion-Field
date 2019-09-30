module Impulse.Impulse where

import           Control.Monad             as M
import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC.CompletionField
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:initDistStr:histFilePath:alphaStr:pinwheelFlagStr:numThreadStr:_) <-
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
      initDist = read initDistStr :: [R2S1RPPoint]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numThread = read numThreadStr :: Int
      sourceDist = L.take 1 initDist
      sinkDist = L.drop 1 initDist
      folderPath = "output/test/Impulse"
      (R2S1RPPoint (_, _, _, s)) = L.head sourceDist
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
      radialArr
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  createDirectoryIfMissing True folderPath
  let arr4d =
        R.slice
          arrR2Z2T0S0
          (Z :. All :. All :. (L.length theta0Freqs - 1) :.
           (L.length scale0Freqs - 1) :.
           All :.
           All)
  createDirectoryIfMissing True (folderPath </> "GreensR2Z2T0S0")
  -- MP.mapM_
  --   (\(i, j) ->
  --      plotImageRepaComplex
  --        (folderPath </> "GreensR2Z2T0S0" </> show (i + 1) L.++ "_" L.++
  --         show (j + 1) L.++
  --         ".png") .
  --      ImageRepa 8 .
  --      computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr4d $
  --      (Z :. i :. j :. All :. All))
  --   [ (i, j)
  --   | i <- [0 .. (L.length thetaFreqs) - 1]
  --   , j <- [0 .. (L.length scaleFreqs) - 1]
  --   ]
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  sourceDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      theta0Freqs
      scale0Freqs
      maxScale
      -- initDist
      sourceDist
  arrR2Z2T0S0F <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ arrR2Z2T0S0
  -- Source field
  sourceArr <- convolveR2T0S0 plan arrR2Z2T0S0F sourceDistArr
  sourceR2Z2 <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    sourceR2Z2
  plotImageRepaComplex (folderPath </> printf "Source%.0f.png" s) . ImageRepa 8 $
    sourceField
