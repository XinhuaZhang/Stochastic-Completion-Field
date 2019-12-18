{-# LANGUAGE FlexibleContexts #-}
module PlotGreensFunctionR2Z2T0S0 where

import           Control.Monad.Parallel  as MP
import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array
import Data.Vector.Unboxed as VU

main = do
  args@(numPointStr:numOrientationStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFileName:alphaStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      thetaFreq = read thetaFreqsStr :: Double
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      theta0Freqs = [-theta0Freq .. -theta0Freq]
      thetaFreqs = [thetaFreq .. thetaFreq]
      scale0Freqs = [(scale0Freq - 1) .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/app/PlotGreensFunctionR2Z2T0S0"
      histPath = folderPath </> histFileName
      radialHistPath = folderPath </> ("Radial_" L.++ histFileName)
  createDirectoryIfMissing True (folderPath </> "GreensR2Z2T0S0")
  flag <- doesFileExist histPath
  hist <-
    if flag
      then decodeFile histPath
      else return $
           emptyHistogram
             [ numPoint
             , numPoint
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0
  arr <-
    solveMonteCarloR2Z2T0S0
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
      histPath
      hist
  let arr4d =
        R.slice
          arr
          (Z :. (L.length theta0Freqs - 1) :. (L.length scale0Freqs - 1) :. All :.
           All :.
           All :.
           All)
      arr4dMax = VU.maximum . VU.map magnitude . toUnboxed . computeS $ arr4d
      -- arr4dNormalized = R.map (\x -> x / (arr4dMax :+ 0) * (255 :+ 0))
  -- MP.mapM_
  --   (\(i, j) ->
  --      plotImageRepaComplex
  --        (folderPath </> "GreensR2Z2T0S0" </> show (i + 1) L.++ "_" L.++
  --         show (j + 1) L.++
  --         ".png") .
  --      ImageRepa 8 .
  --      computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr4dNormalized $
  --      (Z :. i :. j :. All :. All))
  --   [ (i, j)
  --   | i <- [0 .. (L.length thetaFreqs) - 1]
  --   , j <- [0 .. (L.length scaleFreqs) - 1]
  --   ]
  flag <- doesFileExist radialHistPath
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile radialHistPath
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
          radialHistPath
          (emptyHistogram
             [ if maxScale >= 2
                 then round maxScale + 1
                 else (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  arr' <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (PinwheelHollow0 1)
      radialArr
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  let arr4d' =
        R.slice
          arr'
          (Z :. All :. All :. (L.length thetaFreqs - 1) :.
           (L.length scaleFreqs - 1) :.
           All :.
           All)
  MP.mapM_
    (\(i, j) ->
       plotImageRepaComplex
         (folderPath </> "GreensR2Z2T0S0" </> "Radial_" L.++ show (i + 1) L.++
          "_" L.++
          show (j + 1) L.++
          ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr4d' $
       (Z :. i :. j :. All :. All))
    [ (i, j)
    | i <- [0 .. (L.length thetaFreqs) - 1]
    , j <- [0 .. (L.length scaleFreqs) - 1]
    ]
