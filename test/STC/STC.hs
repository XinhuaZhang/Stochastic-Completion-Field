{-# LANGUAGE BangPatterns #-}
module STC where

import           Control.Monad.Parallel      as MP
import           Data.Array.IArray           as IA
import           Data.Array.Repa             as R
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.GreensFunction
import           FokkerPlanck.Histogram
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Parallel
import           Utils.Time

main = do
  args@(gpuIDStr:numPointStr:deltaStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:deltaLogStr:taoStr:numTrailStr:maxTrailStr:r2FreqStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initDistStr:initScaleStr:histFilePath:stdStr:numThreadStr:_) <-
    getArgs
  let gpuID = read gpuIDStr :: [Int]
      numPoint = read numPointStr :: Int
      delta = read deltaStr :: Double
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      tao = read taoStr :: Double
      r2Freq = read r2FreqStr :: Int
      phiFreq = read phiFreqsStr :: Int
      phiFreqs = L.map fromIntegral [-phiFreq .. phiFreq]
      rhoFreq = read rhoFreqsStr :: Int
      rhoFreqs = L.map fromIntegral [-rhoFreq .. rhoFreq]
      thetaFreq = read thetaFreqsStr :: Int
      thetaFreqs = L.map fromIntegral [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Int
      scaleFreqs = L.map fromIntegral [-scaleFreq .. scaleFreq]
      initScale = read initScaleStr :: Double
      initDist = read initDistStr :: [(Double, Double, Double, Double)]
      initPoints = L.map (\(x, y, t, s) -> Point x y t s) initDist
      numThread = read numThreadStr :: Int
      folderPath = "output/test/FourierSeriesOfPinwheels"
      maxScale = read maxScaleStr :: Double
      halfLogPeriod = log maxScale
      deltaLog = read deltaLogStr :: Double
      std = read stdStr :: Double
      radius = fromIntegral $ div numPoint 2
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  printCurrentTime "Start computing coefficients..."
  hist <-
    sampleCartesian
      folderPath
      histFilePath
      gpuID
      radius
      144
      0.25 -- deltaLog
      initScale
      thetaSigma
      tao
      phiFreqs
      rhoFreqs
      thetaFreqs
      scaleFreqs
  printCurrentTime "Done"
  printCurrentTime "Start Convloution.."
  let !coefficients =
        normalizeFreqArr' std phiFreqs rhoFreqs . getNormalizedHistogramArr $
        hist
      initialDistribution =
        computeInitialDistributionFull' r2Freq radius phiFreq rhoFreq initPoints
      !harmonicsArray' =
        computeFourierSeriesOfLogPolarHarmonicsArray
          radius
          deltaLog
          r2Freq
          phiFreq
          rhoFreq
          thetaFreq
          scaleFreq
          halfLogPeriod :: IA.Array (Int, Int) (VS.Vector (Complex Double))
      -- !harmonicsArray' = computeHarmonicsArray numPoint 1 numPoint 1 phiFreqs rhoFreqs thetaFreqs scaleFreqs 1 radius
      source = convolveFull' coefficients harmonicsArray' initialDistribution
      inverseR2Harmonics =
        computeRectangularInverseHarmonics numPoint numPoint delta radius r2Freq :: [VS.Vector (Complex Double)]
      -- pinwheelVec =
      --   parMap
      --     rdeepseq
      --     (\(idx, vec) ->
      --        ( idx
      --        , VU.fromList . L.map (VS.sum . VS.zipWith (*) vec) $
      --          inverseR2Harmonics)) .
      --   assocs $
      --   harmonicsArray'
  -- print initialDistribution
  createDirectoryIfMissing True (folderPath </> "Pinwheels")
  plotImageRepaComplex (folderPath </> "Pinwheels" </> "Pinwheel_(1,1).png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) .
    VU.fromList . L.map (VS.sum . VS.zipWith (*) (harmonicsArray' IA.! (1, 1))) $
    inverseR2Harmonics
  -- MP.mapM_
  --   (\((a, b), vec) ->
  --      plotImageRepaComplex
  --        (folderPath </> "Pinwheels" </> printf "Pinwheel_(%d,%d).png" a b) .
  --      ImageRepa 8 . fromUnboxed (Z :. (1 :: Int) :. numPoint :. numPoint) $
  --      vec)
  --   pinwheelVec
  -- printCurrentTime "Done"
  printCurrentTime "Start ploting.."
  plotDFTArrayPower
    (folderPath </> "SourcePower.png")
    numPoint
    numPoint
    -- (2 * r2Freq + 1)
    -- (2 * r2Freq + 1)
    source
  plotFourierCoefficientsPower
    (folderPath </> "Source.png")
    numPoint
    numPoint
    inverseR2Harmonics
    source
  printCurrentTime "Done"
