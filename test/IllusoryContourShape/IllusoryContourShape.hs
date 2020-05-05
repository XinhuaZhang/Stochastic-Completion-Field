{-# LANGUAGE BangPatterns #-} 
module IllusoryContourShape where

import           Control.Monad               as M
import           Data.Array.Repa             as R
import           Data.Binary
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Generic         as VG
import           FokkerPlanck
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.Time
import qualified Data.Array.IArray as IA
import           FokkerPlanck.GreensFunction

main = do
  args@(gpuIDStr:numPointStr:deltaXStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:cutoffStr:deltaLogStr:taoStr:numTrailStr:maxTrailStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:initScaleStr:shape2DStr:histFilePath:writeFlagStr:numIterationStr:suffix:stdStr:stdGStr:numThreadStr:_) <-
    getArgs
  let gpuID = read gpuIDStr :: [Int]
      numPoint = read numPointStr :: Int
      deltaX = read deltaXStr :: Double
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      phiFreq = read phiFreqsStr :: Double
      phiFreqs = [-phiFreq .. phiFreq]
      rhoFreq = read rhoFreqsStr :: Double
      rhoFreqs = [-rhoFreq .. rhoFreq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Double
      scaleFreqs = [-scaleFreq .. scaleFreq]
      initScale = read initScaleStr :: Double
      shape2D@(Points _ minDist _) = read shape2DStr :: Points Shape2D
      numThread = read numThreadStr :: Int
      folderPath =
        "output/test/IllusoryContourShape" </> (takeBaseName histFilePath) L.++
        "_" L.++
        cutoffStr L.++
        "_" L.++
        deltaXStr
      maxScale = read maxScaleStr :: Double
      cutoff = read cutoffStr :: Double
      halfLogPeriod = log maxScale
      deltaLog = read deltaLogStr :: Double
      writeFlag = read writeFlagStr :: Bool
      std = read stdStr :: Double
      stdG = read stdGStr :: Double
      numIteration = read numIterationStr :: Int
  print thetaFreqs
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  hist <-
    if flag
      then do
        printCurrentTime $
          "read Fourier coefficients data from " L.++ histFilePath
        decodeFile histFilePath
      -- else runMonteCarloFourierCoefficientsGPU
      --        gpuID
      --        numThread
      --        numTrail
      --        maxTrail
      --        thetaSigma
      --        scaleSigma
      --        maxScale
      --        tao
      --        phiFreqs
      --        rhoFreqs
      --        thetaFreqs
      --        scaleFreqs
      --        deltaLog
      --        initScale
      --        histFilePath
      --        (emptyHistogram
      --           [ L.length phiFreqs
      --           , L.length rhoFreqs
      --           , L.length thetaFreqs
      --           , L.length scaleFreqs
      --           ]
      --           0)
      else sampleCartesian
             folderPath
             histFilePath
             gpuID
             numPoint
             numPoint
             180
             deltaX
             deltaLog
             initScale
             thetaSigma
             tao
             phiFreqs
             rhoFreqs
             thetaFreqs
             scaleFreqs
  let !points =
        L.map (\(!x, !y) -> Point x y 0 1) . getShape2DIndexList . makeShape2D $
        shape2D
      -- !xs = L.map (\(a, b) -> (a * deltaX, b * deltaX)) . makeShape2D $ shape2D
      !xs = getShape2DIndexList . makeShape2D $ shape2D
      !coefficients =
        normalizeFreqArr std phiFreqs rhoFreqs . getNormalizedHistogramArr $
        hist
       -- = getNormalizedHistogramArr $ hist :: R.Array U DIM4 (Complex Double)
      !thetaRHarmonics =
        computeThetaRHarmonics
          numOrientation
          numScale
          thetaFreqs
          scaleFreqs
          halfLogPeriod
  print xs
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPoint
      numPoint
      (L.length thetaFreqs)
      (L.length scaleFreqs)
  let !initSourceSparse =
        computeInitialDistributionPowerMethodSparse phiFreqs rhoFreqs points
      !harmonicsArraySparse =
        computeHarmonicsArraySparse
          numPoint
          deltaX
          numPoint
          deltaX
          phiFreqs
          rhoFreqs
          thetaFreqs
          scaleFreqs
          halfLogPeriod
          cutoff
  -- let !initSource =
  --       computeInitialDistributionPowerMethod'
  --         numPoint
  --         numPoint
  --         phiFreqs
  --         rhoFreqs
  --         -- thetaFreqs
  --         -- scaleFreqs
  --         points
  --     !bias = computeBias numPoint numPoint points
  gaussian <- gaussianFilter2D plan numPoint stdG
  -- harmonicsArray <-
  --   dftHarmonicsArrayG
  --     plan
  --     numPoint
  --     deltaX
  --     numPoint
  --     deltaX
  --     phiFreqs
  --     rhoFreqs
  --     thetaFreqs
  --     scaleFreqs
  --     halfLogPeriod
  --     cutoff
  --     gaussian
  -- print . IA.indices $ harmonicsArraySparse
  -- completion <-
  computeContourSparse
    plan
    folderPath
    coefficients
    harmonicsArraySparse
    -- thetaRHarmonics
    (fromListUnboxed (Z :. L.length phiFreqs) phiFreqs)
    (fromListUnboxed (Z :. L.length rhoFreqs) rhoFreqs)
    cutoff
    -- gaussian
                      -- std
    xs
    numIteration
    suffix
    initSourceSparse
  -- completion <-
  --   computeContour'
  --    plan
  --    folderPath
  --    writeFlag
  --    coefficients
  --    harmonicsArray
  --    bias
  --    -- gaussian
  --    numIteration
  --    suffix
  --    thetaRHarmonics
  --    initSource
  -- plotDFTArrayThetaR
  --   (folderPath </> (printf "Completion_%s.png" suffix))
  --   numPoint
  --   numPoint
  --   thetaRHarmonics
  --   completion
  -- plotDFTArrayThetaRMag
  --   (folderPath </> (printf "CompletionMax_%s.png" suffix))
  --   numPoint
  --   numPoint
  --   thetaRHarmonics
  --   completion
  -- computeContourSparse'''
  --   plan
  --   folderPath
  --   coefficients
  --   harmonicsArraySparse
  --   thetaRHarmonics
  --   (fromListUnboxed (Z :. L.length phiFreqs) phiFreqs)
  --   (fromListUnboxed (Z :. L.length rhoFreqs) rhoFreqs)
  --   cutoff
  --   xs
  --   numIteration
  --   (suffix L.++ "Test")
  --   initSourceSparse
