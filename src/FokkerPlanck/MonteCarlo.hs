{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.MonteCarlo
  ( runMonteCarloFourierCoefficients
  , runMonteCarloFourierCoefficientsGPU
  , solveMonteCarloR2S1
  ) where

import           Array.UnboxedArray              as UA
import           Control.DeepSeq
import           Control.Monad                   as M
import           Control.Monad.Parallel          as MP
import           Data.Array.Accelerate           (constant)
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                 as R
import           Data.Binary
import           Data.ByteString.Lazy            as BL
import           Data.Complex
import           Data.DList                      as DL
import           Data.Ix
import           Data.List                       as L
import           Data.Vector.Unboxed             as VU
import           FokkerPlanck.BrownianMotion
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.FourierSeriesGPU
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver             as CUDA
import           Statistics.Distribution.Normal
import           Statistics.Distribution.Uniform
import           System.Directory
import           System.IO                       as IO
import           System.Random.MWC
import           Text.Printf
import           Utils.Parallel
import           Utils.Time

{-# INLINE runBatch #-}
runBatch :: Int -> Int -> Int -> (Int -> a -> [GenIO] -> IO a) -> a -> IO a
runBatch !numGens !numTrails !batchSize monterCarloHistFunc !init = do
  let (!numBatch, !numLeftover) = divMod numTrails batchSize
  printf
    "%d trails and batch size is %d, %d batch in total.\n"
    numTrails
    batchSize
    (numBatch +
     if numLeftover > 0
       then 1
       else 0)
  gensList <- M.replicateM numBatch (M.replicateM numGens createSystemRandom)
  x <-
    M.foldM
      (\y gens -> do
         printCurrentTime ""
         monterCarloHistFunc batchSize y gens)
      init
      gensList
  if numLeftover > 0
    then do
      gens <- M.replicateM numGens createSystemRandom
      monterCarloHistFunc numLeftover x gens
    else return x

{-# INLINE computeHistogramFromMonteCarloParallel #-}
computeHistogramFromMonteCarloParallel ::
     (NFData hist, Binary hist, Show hist)
  => FilePath
  -> (GenIO -> IO particle)
  -> ([particle] -> hist)
  -> (hist -> hist -> hist)
  -> Int
  -> hist
  -> [GenIO]
  -> IO hist
computeHistogramFromMonteCarloParallel !filePath pointsGenerator histFunc addHist !n !initHist !gens = do
  xs <- MP.mapM (M.replicateM (div n . L.length $ gens) . pointsGenerator) gens
  let tmpFilePath = (filePath L.++ "_tmp")
      hist = L.foldl' addHist initHist . parMap rdeepseq histFunc $ xs
  unless
    (L.null filePath)
    (do encodeFile tmpFilePath hist
        copyFile tmpFilePath filePath)
  return hist

{-# INLINE computeHistogramFromMonteCarloParallelSingleGPU #-}
computeHistogramFromMonteCarloParallelSingleGPU ::
     (NFData hist, Binary hist, Show hist)
  => FilePath
  -> (GenIO -> IO particle)
  -> ([particle] -> hist)
  -> (hist -> hist -> hist)
  -> Int
  -> hist
  -> [GenIO]
  -> IO hist
computeHistogramFromMonteCarloParallelSingleGPU !filePath pointsGenerator histFunc addHist !n !initHist !gens = do
  xs <- MP.mapM (M.replicateM (div n . L.length $ gens) . pointsGenerator) gens
  let tmpFilePath = (filePath L.++ "_tmp")
      hist = addHist initHist . histFunc . L.concat $ xs
  unless
    (L.null filePath)
    (do encodeFile tmpFilePath hist
        copyFile tmpFilePath filePath)
  return hist
  
{-# INLINE computeHistogramFromMonteCarloParallelMultipleGPU #-}
computeHistogramFromMonteCarloParallelMultipleGPU ::
     (NFData hist, Binary hist, Show hist)
  => FilePath
  -> (GenIO -> IO particle)
  -> (PTX -> [particle] -> hist)
  -> (hist -> hist -> hist)
  -> [PTX]
  -> Int
  -> hist
  -> [GenIO]
  -> IO hist
computeHistogramFromMonteCarloParallelMultipleGPU !filePath pointsGenerator histFunc addHist !ctxs !n !initHist !gens = do
  xs <- MP.mapM (M.replicateM (div n . L.length $ gens) . pointsGenerator) gens
  let tmpFilePath = (filePath L.++ "_tmp")
      hist = L.foldl' addHist initHist . parZipWith rdeepseq histFunc ctxs $ xs
  unless
    (L.null filePath)
    (do encodeFile tmpFilePath hist
        copyFile tmpFilePath filePath)
  return hist

{-# INLINE runMonteCarloFourierCoefficients #-}
runMonteCarloFourierCoefficients ::
     Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (Histogram (Complex Double))
runMonteCarloFourierCoefficients !numGens !numTrails !batchSize !thetaSigma !scaleSigma !maxScale !tao !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !deltaLog !filePath !initHist = do
  when
    (maxScale <= 1)
    (error $
     printf "runMonteCarloFourierCoefficients error: maxScale(%f) <= 1" maxScale)
  let !thetaDist = normalDistrE 0 thetaSigma
      !scaleDist = normalDistrE 0 scaleSigma
      pointsGenerator = generatePath thetaDist scaleDist maxScale tao
      histFunc =
        computeFourierCoefficients
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          (log maxScale)
          deltaLog
      monterCarloHistFunc =
        computeHistogramFromMonteCarloParallel
          filePath
          pointsGenerator
          histFunc
          addHistogramUnsafe
  hist <- runBatch numGens numTrails batchSize monterCarloHistFunc initHist
  return hist


{-# INLINE runMonteCarloFourierCoefficientsGPU #-}
runMonteCarloFourierCoefficientsGPU ::
     [Int]
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (Histogram (Complex Double))
runMonteCarloFourierCoefficientsGPU !deviceIDs !numGens !numTrails !batchSize !thetaSigma !scaleSigma !maxScale !tao !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !deltaLog !filePath !initHist = do
  when
    (maxScale <= 1)
    (error $
     printf "runMonteCarloFourierCoefficients error: maxScale(%f) <= 1" maxScale)
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let !thetaDist = normalDistrE 0 thetaSigma
      !scaleDist = normalDistrE 0 scaleSigma
      !freqArr = computeFrequencyArray phiFreqs rhoFreqs thetaFreqs rFreqs
      pointsGenerator = generatePath thetaDist scaleDist maxScale tao
      histFuncSingleGPU =
        computeFourierCoefficientsGPU
          freqArr
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          (log maxScale)
          deltaLog
          (constant maxScale)
          (constant (log maxScale))
          (constant (deltaLog :+ 0))
          (L.head ptxs)
      histFuncMultipleGPU =
        computeFourierCoefficientsGPU
          freqArr
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          (log maxScale)
          deltaLog
          (constant maxScale)
          (constant (log maxScale))
          (constant (deltaLog :+ 0))
      monterCarloHistFuncSingleGPU =
        computeHistogramFromMonteCarloParallelSingleGPU
          filePath
          pointsGenerator
          histFuncSingleGPU
          addHistogramUnsafe
      monterCarloHistFuncMultipleGPU =
        computeHistogramFromMonteCarloParallelMultipleGPU
          filePath
          pointsGenerator
          histFuncMultipleGPU
          addHistogramUnsafe
          ptxs
  hist <-
    runBatch
      numGens
      numTrails
      batchSize
      (if L.length deviceIDs == 1
         then monterCarloHistFuncSingleGPU
         else monterCarloHistFuncMultipleGPU)
      initHist
  M.mapM destroy ctxs
  return hist

{-# INLINE countR2S1 #-}
countR2S1 :: (Int, Int) -> (Int, Int) -> Int -> [DList Particle] -> Histogram Double
countR2S1 xRange@(!xMin, !xMax) yRange@(!yMin, !yMax) !numOrientations !xs =
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        L.filter (\((_, x, y), _) -> (inRange xRange x) && (inRange yRange y)) .
        L.map
          (\(Particle phi rho theta _) ->
             let !x = rho * cos phi
                 !y = rho * sin phi
             in ((floor $ theta / deltaTheta, round x, round y), 1)) .
        DL.toList . DL.concat $
        xs
      !numTrajectories = L.length ys
      !arr =
        UA.accum (+) 0 ((0, xMin, yMin), (numOrientations - 1, xMax, yMax)) ys
  in Histogram
       [(yMax - yMin + 1), (xMax - xMin + 1), numOrientations]
       numTrajectories
       (toUnboxedVector arr)

{-# INLINE solveMonteCarloR2S1 #-}
solveMonteCarloR2S1 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> FilePath
  -> IO (R.Array U DIM3 Double)
solveMonteCarloR2S1 numGens numTrails xLen yLen numOrientations thetaSigma tao r histFilePath = do
  let !xShift = div xLen 2
      xRange' =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      yRange' =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 0
      pointsGenerator =
        generatePath
          thetaDist
          scaleDist
          (sqrt . fromIntegral $ xLen ^ 2 + yLen ^ 2)
          tao
      histFunc = countR2S1 xRange' yRange' numOrientations
      monterCarloHistFunc =
        computeHistogramFromMonteCarloParallel
          ""
          pointsGenerator
          histFunc
          addHistogramUnsafe
  hist <-
    runBatch
      numGens
      numTrails
      numTrails
      monterCarloHistFunc
      (emptyHistogram [yLen, xLen, numOrientations] 0)
  return . getNormalizedHistogramArr $ hist
