{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.MonteCarlo
  ( -- runMonteCarloFourierCoefficients
    runMonteCarloFourierCoefficientsGPU
  , solveMonteCarloR2S1
  ) where

import           Array.UnboxedArray              as UA
import           Control.DeepSeq
import           Control.Monad                   as M
import           Control.Concurrent.Async
import           Control.Monad.Parallel          as MP
import           Data.Array.Accelerate           (constant, use)
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
import           Statistics.Distribution
import           Statistics.Distribution.Laplace
import           Statistics.Distribution.Normal
import           Statistics.Distribution.Poisson
import           Statistics.Distribution.Uniform
import           System.Directory
import           System.IO                       as IO
import           System.Random.MWC
import           Text.Printf
import           Utils.List
import           Utils.Parallel
import           Utils.Time
import GHC.Float

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
  -> ([particle] -> IO hist)
  -> (hist -> hist -> hist)
  -> Int
  -> hist
  -> [GenIO]
  -> IO hist
computeHistogramFromMonteCarloParallel !filePath pointsGenerator histFunc addHist !n !initHist !gens = do
  xs <- MP.mapM (M.replicateM (div n . L.length $ gens) . pointsGenerator) gens
  let tmpFilePath = (filePath L.++ "_tmp")
  hist <- L.foldl' addHist initHist <$> M.mapM histFunc xs
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
  -> Int
  -> hist
  -> [GenIO]
  -> IO hist
computeHistogramFromMonteCarloParallelSingleGPU !filePath pointsGenerator histFunc addHist !deviceID !n !initHist !gens = do
  xs <- MP.mapM (M.replicateM (div n . L.length $ gens) . pointsGenerator) gens
  let tmpFilePath = (filePath L.++ "_tmp")
      !hist = addHist initHist . histFunc . L.concat $ xs
  unless
    (L.null filePath)
    (do encodeFile tmpFilePath hist
        copyFile tmpFilePath filePath)
  return hist

{-# INLINE computeHistogramFromMonteCarloParallelMultipleGPU #-}
computeHistogramFromMonteCarloParallelMultipleGPU ::
     (NFData hist, Binary hist, Show hist, NFData particle)
  => FilePath
  -> (GenIO -> IO particle)
  -> (PTX -> [particle] -> hist)
  -> (hist -> hist -> hist)
  -> [PTX]
  -> Int
  -> hist
  -> [GenIO]
  -> IO hist
computeHistogramFromMonteCarloParallelMultipleGPU !filePath pointsGenerator histFunc addHist !ptxs !n !initHist !gens = do
  xs <-
    L.concat <$>
    mapConcurrently
      (M.replicateM (div n . L.length $ gens) . pointsGenerator)
      gens
  let tmpFilePath = filePath L.++ "_tmp"
      !hist =
        L.foldl' addHist initHist .
        parZipWith rdeepseq histFunc ptxs . divideListN (L.length ptxs) $
        xs
  unless
    (L.null filePath)
    (do encodeFile tmpFilePath hist
        copyFile tmpFilePath filePath)
  return hist

-- {-# INLINE runMonteCarloFourierCoefficients #-}
-- runMonteCarloFourierCoefficients ::
--      Int
--   -> Int
--   -> Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> FilePath
--   -> Histogram (Complex Double)
--   -> IO (Histogram (Complex Double))
-- runMonteCarloFourierCoefficients !numGens !numTrails !batchSize !thetaSigma !scaleSigma !maxScale !tao !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !deltaLog !filePath !initHist = do
--   when
--     (maxScale <= 1)
--     (error $
--      printf "runMonteCarloFourierCoefficients error: maxScale(%f) <= 1" maxScale)
--   let !thetaDist = normalDistrE 0 thetaSigma
--       !scaleDist = normalDistrE 0 scaleSigma
--       pointsGenerator = generatePath thetaDist scaleDist maxScale tao 1
--       histFunc =
--         computeFourierCoefficients
--           phiFreqs
--           rhoFreqs
--           thetaFreqs
--           rFreqs
--           (log maxScale)
--           deltaLog
--       monterCarloHistFunc =
--         computeHistogramFromMonteCarloParallel
--           filePath
--           pointsGenerator
--           histFunc
--           addHistogramUnsafe
--   hist <- runBatch numGens numTrails batchSize monterCarloHistFunc initHist
--   return hist


{-# INLINE runMonteCarloFourierCoefficientsSingleGPU #-}
runMonteCarloFourierCoefficientsSingleGPU ::
     Int
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
  -> Double -> Double
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (Histogram (Complex Double))
runMonteCarloFourierCoefficientsSingleGPU !deviceID !numGens !numTrails !batchSize !thetaSigma !scaleSigma !maxScale !tao !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !deltaLog !initScale !filePath !initHist = do
  when
    (maxScale <= 1)
    (error $
     printf "runMonteCarloFourierCoefficients error: maxScale(%f) <= 1" maxScale)
  gen <- createSystemRandom
  initialise []
  dev <- device deviceID
  ctx <- CUDA.create dev []
  ptx <- createTargetFromContext ctx
  let !thetaDist = normalDistrE 0 thetaSigma
      !scaleDist = normalDistrE 0 scaleSigma
      -- !scaleDist = uniformDistrE (- (log maxScale)) (log maxScale)
      !poissonDist = poissonE 0
      !freqArr =
        computeFrequencyArray
          (L.map double2Float phiFreqs)
          (L.map double2Float rhoFreqs)
          (L.map double2Float thetaFreqs)
          (L.map double2Float rFreqs)
      pointsGenerator =
        generatePath
          thetaDist
          scaleDist
          poissonDist
          (maxScale ^ 2)
          maxScale
          tao
          initScale
          1
      histFuncSingleGPU =
        computeFourierCoefficientsGPU
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          (use freqArr)
          (constant . double2Float $ maxScale)
          (constant . double2Float . log $ maxScale)
          ptx
      monterCarloHistFuncSingleGPU =
        computeHistogramFromMonteCarloParallelSingleGPU
          filePath
          pointsGenerator
          histFuncSingleGPU
          addHistogramUnsafe
          deviceID
  hist <-
    runBatch numGens numTrails batchSize monterCarloHistFuncSingleGPU initHist
  -- deepseq hist (destroy ctx)
  return hist

-- {-# INLINE runMonteCarloFourierCoefficientsMultipleGPU #-}
-- runMonteCarloFourierCoefficientsMultipleGPU ::
--      [Int]
--   -> Int
--   -> Int
--   -> Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> Double
--   -> FilePath
--   -> Histogram (Complex Double)
--   -> IO (Histogram (Complex Double))
-- runMonteCarloFourierCoefficientsMultipleGPU !deviceIDs !numGens !numTrails !batchSize !thetaSigma !scaleSigma !maxScale !tao !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !deltaLog !initScale !filePath !initHist = do
--   when
--     (maxScale <= 1)
--     (error $
--      printf "runMonteCarloFourierCoefficients error: maxScale(%f) <= 1" maxScale)
--   initialise []
--   devs <- M.mapM device deviceIDs
--   ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
--   ptxs <- M.mapM createTargetFromContext ctxs
--   let !thetaDist = normalDistrE 0 thetaSigma
--       !scaleDist = normalDistrE 0 scaleSigma
--       !freqArr = computeFrequencyArray phiFreqs rhoFreqs thetaFreqs rFreqs
--       pointsGenerator = generatePath thetaDist scaleDist maxScale tao initScale
--       histFuncMultipleGPU =
--         computeFourierCoefficientsGPU
--           phiFreqs
--           rhoFreqs
--           thetaFreqs
--           rFreqs
--           (log maxScale)
--           deltaLog
--           (use freqArr)
--           (constant maxScale)
--           (constant (log maxScale))
--           (constant (deltaLog :+ 0))
--       monterCarloHistFuncMultipleGPU =
--         computeHistogramFromMonteCarloParallelMultipleGPU
--           filePath
--           pointsGenerator
--           histFuncMultipleGPU
--           addHistogramUnsafe
--           ptxs
--   hist <-
--     runBatch numGens numTrails batchSize monterCarloHistFuncMultipleGPU initHist
--   deepseq hist (M.mapM destroy ctxs)
--   return hist

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
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> FilePath
  -> IO (Histogram (Complex Double))
runMonteCarloFourierCoefficientsGPU !deviceIDs !numGens !numTrails !batchSize !thetaSigma !scaleLambda !poissonLambda !sigma !tao !deltaT !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !periodEnv !stdR2 !filePath = do
  when
    (periodEnv <= 1)
    (error $
     printf
       "runMonteCarloFourierCoefficients error: periodEnv(%f) <= 1"
       periodEnv)
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let !initHist =
        emptyHistogram
          [ L.length phiFreqs
          , L.length rhoFreqs
          , L.length thetaFreqs
          , L.length rFreqs
          ]
          0
      !thetaDist = normalDistrE 0 (thetaSigma * sqrt deltaT)
      !scaleDist = normalDistrE 0 (scaleLambda * sqrt deltaT)
      !poissonDist = poissonE (poissonLambda / deltaT)
      !freqArr =
        computeFrequencyArray
          (L.map double2Float phiFreqs)
          (L.map double2Float rhoFreqs)
          (L.map double2Float thetaFreqs)
          (L.map double2Float rFreqs)
      !maxRho =  sqrt periodEnv / 1 / sqrt 2
      !maxR = sqrt periodEnv / 1 / sqrt 2
      pointsGenerator =
        generatePath
          thetaDist
          scaleDist
          poissonDist
          maxRho
          maxR
          (tao / deltaT)
          deltaT
          stdR2
      -- pointsGenerator =
      --   generatePath'
      --     thetaSigma
      --     scaleDist
      --     poissonLambda
      --     maxRho
      --     maxR
      --     tao
      --     deltaT
      --     stdR2
      -- pointsGenerator =
      --   generatePath' thetaSigma scaleDist poissonLambda periodEnv (tao / deltaT) deltaT
      histFuncSingleGPU =
        computeFourierCoefficientsGPU
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          (use freqArr)
          (constant . double2Float $ sigma)
          (constant . double2Float $ log periodEnv)
      monterCarloHistFuncGPU =
        computeHistogramFromMonteCarloParallelMultipleGPU
          filePath
          pointsGenerator
          histFuncSingleGPU
          addHistogramUnsafe
          ptxs
  runBatch numGens numTrails batchSize monterCarloHistFuncGPU initHist

{-# INLINE countR2S1 #-}
countR2S1 ::
     GenIO
  -> Double
  -> (Int, Int)
  -> (Int, Int)
  -> Int
  -> [DList Particle]
  -> IO (Histogram Double)
countR2S1 randomGen thetaSigma xRange@(!xMin, !xMax) yRange@(!yMin, !yMax) !numOrientations !xs = do
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
  ys <-
    fmap
      (L.filter
         (\((_, x, y), _) ->
            (inRange xRange x) && (inRange yRange y) && (x /= 0 || y /= 0)) .
       L.map
         (\(Particle phi rho theta' _ _) ->
            let !x = rho * cos phi
                !y = rho * sin phi
                !theta =
                  if theta' < 0
                    then theta' + 2 * pi
                    else theta'
             in ((floor $ theta / deltaTheta, round x, round y), 1)) .
       L.concat) .
    M.mapM
      (\particle@(Particle phi rho theta r _) -> do
         let delta = 1
             n = Prelude.floor $ r / delta
         zs <-
           M.mapM
             (\i -> do
                let thetaDist =
                      normalDistr
                        0
                        (thetaSigma * sqrt (delta * fromIntegral i / r))
                    (Particle a b c d _) =
                      FokkerPlanck.BrownianMotion.moveParticle
                        1
                        (Particle
                           phi
                           rho
                           theta
                           (Prelude.fromIntegral i * delta)
                           1)
                deltaThetaDiffusion <- genContVar thetaDist randomGen
                return (Particle a b (c `thetaPlus` deltaThetaDiffusion) d 1))
             [1 .. n - 1]
         return .
           L.concatMap
             (\(Particle phi' rho' theta' r' v) ->
                [ (Particle phi' rho' theta' r' v)
                , (Particle (-phi') rho' (-theta') r' v)
                ]) $
           (FokkerPlanck.BrownianMotion.moveParticle
              1
              (Particle phi rho theta r 1)) :
           zs) .
    DL.toList . DL.concat $
    xs
  return .
    Histogram
      [(yMax - yMin + 1), (xMax - xMin + 1), numOrientations]
      (L.length ys) .
    toUnboxedVector .
    UA.accum (+) 0 ((0, xMin, yMin), (numOrientations - 1, xMax, yMax)) $
    ys

{-# INLINE solveMonteCarloR2S1 #-}
solveMonteCarloR2S1 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> FilePath
  -> IO (R.Array U DIM3 Double)
solveMonteCarloR2S1 numGens numTrails batchSize xLen yLen numOrientations thetaSigma tao r initSpeed histFilePath = do
  gen <- createSystemRandom
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
      poissonDist = poissonE 0
      pointsGenerator =
        generatePath
          thetaDist
          scaleDist
          poissonDist
          (r^2)
          r -- (sqrt . fromIntegral $ xLen ^ 2 + yLen ^ 2)
          tao
          initSpeed
          1
      histFunc = countR2S1 gen thetaSigma xRange' yRange' numOrientations
      monterCarloHistFunc =
        computeHistogramFromMonteCarloParallel
          histFilePath
          pointsGenerator
          histFunc
          addHistogramUnsafe
  hist <-
    runBatch
      numGens
      numTrails
      batchSize
      monterCarloHistFunc
      (emptyHistogram [yLen, xLen, numOrientations] 0)
  return . getNormalizedHistogramArr $ hist
