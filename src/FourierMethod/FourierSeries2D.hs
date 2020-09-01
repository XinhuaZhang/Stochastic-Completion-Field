{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
module FourierMethod.FourierSeries2D
  ( module FourierMethod.FourierSeries2D
  , module FourierMethod.BlockCudaMatrix
  , module FourierMethod.BlockMatrixAcc
  ) where

import           Control.DeepSeq
import           Control.Monad                               as M
import           Control.Monad.IO.Class
import           Control.Monad.Trans.Resource
import qualified Data.Array.Accelerate                       as A
import           Data.Array.Accelerate.LLVM.PTX              as A
import           Data.Array.Accelerate.Numeric.LinearAlgebra as A
import           Data.Array.Repa                             as R
import           Data.Complex
import           Data.Conduit                                as C
import           Data.Conduit.List                           as CL
import           Data.List                                   as L
import           Data.Vector.Storable                        as VS
import           Data.Vector.Unboxed                         as VU
import           Foreign.CUDA.Driver                         as CUDA
import           FourierMethod.BlockCudaMatrix
import           FourierMethod.BlockMatrixAcc
import           FourierMethod.FourierSeries2DAcc
import           Utils.BLAS
import           Utils.Distribution
import           Utils.List
import           Utils.Parallel                              hiding ((.|))
import           Utils.Time

{-# INLINE harmonicMatPTX #-}
harmonicMatPTX ::
     ( Storable e
     , A.Elt (Complex e)
     , A.Elt e
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => Int
  -> e
  -> e
  -> e
  -> PTX
  -> [(Int, Int)]
  -> CuMat (Complex e)
harmonicMatPTX numPoints period delta deltaFreq ptx r2Freqs =
  let rows = L.length r2Freqs
  in CuMat rows (numPoints ^ 2) .
     CuVecHost .
     VS.fromList .
     A.toList .
     run1With ptx (harmonicAcc numPoints period delta deltaFreq) .
     A.fromList (A.Z A.:. (L.length r2Freqs)) $
     r2Freqs

-- create column-major inverse harmonics
{-# INLINE inverseHarmonicMatPTX #-}
inverseHarmonicMatPTX ::
     ( Storable e
     , A.Elt (Complex e)
     , A.Elt e
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => Int
  -> e
  -> e
  -> PTX
  -> [(Int, Int)]
  -> CuMat (Complex e)
inverseHarmonicMatPTX numFreqs period delta ptx r2Positions =
  let cols = L.length r2Positions
  in CuMat (numFreqs ^ 2) cols .
     CuVecHost .
     VS.fromList .
     A.toList .
     run1With ptx (inverseHarmonicAcc numFreqs period delta) .
     A.fromList (A.Z A.:. (L.length r2Positions)) $
     r2Positions

{-# INLINE inverseHarmonicMatPTX1 #-}
inverseHarmonicMatPTX1 ::
     ( Storable e
     , A.Elt (Complex e)
     , A.Elt e
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => Int
  -> e
  -> e
  -> PTX
  -> [(Int, Int)]
  -> CuMat (Complex e)
inverseHarmonicMatPTX1 numFreqs period delta ptx r2Positions =
  let rows = L.length r2Positions
  in CuMat rows (numFreqs ^ 2) .
     CuVecHost .
     VS.fromList .
     A.toList .
     run1With ptx (inverseHarmonicAcc1 numFreqs period delta) .
     A.fromList (A.Z A.:. (L.length r2Positions)) $
     r2Positions

createHarmonicMatriesGPU ::
     ( Storable e
     , A.Elt (Complex e)
     , A.Elt e
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => [PTX]
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> IO [[CuMat (Complex e)]]
createHarmonicMatriesGPU ptxs numBatch numPoints numFreqs period delta deltaFreq = do
  let idxs =
        L.map (divideListN numBatch) . divideListN (L.length ptxs) $
        [ (xFreq, yFreq)
        | xFreq <- getListFromNumber numFreqs
        , yFreq <- getListFromNumber numFreqs
        ]
  let output =
        parZipWith
          rdeepseq
          (\ptx -> L.map (harmonicMatPTX numPoints period delta deltaFreq ptx))
          ptxs $
        idxs
  return output


createInverseHarmonicMatriesGPU ::
     ( Storable e
     , A.Elt (Complex e)
     , A.Elt e
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => [PTX]
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> IO [[CuMat (Complex e)]]
createInverseHarmonicMatriesGPU ptxs numBatch numPoints numFreqs period delta = do
  let idxs =
        L.map (divideListN numBatch) . divideListN (L.length ptxs) $
        [ (x, y)
        | x <- getListFromNumber numPoints
        , y <- getListFromNumber numPoints
        ]
  return .
    parZipWith
      rdeepseq
      (\ptx -> L.map (inverseHarmonicMatPTX numFreqs period delta ptx))
      ptxs $
    idxs

computeFourierCoefficientsR2 ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> Int
  -> [CuMat (Complex e)]
  -> IO (R.Array U DIM3 (Complex e))
computeFourierCoefficientsR2 deviceIDs ptxs numFreqs numPoints period delta deltaFreq numBatchR2Freqs xs = do
  let simpsonNorm = (delta / 3) ^ 2 :+ 0
      cols = L.foldl' (\s mat -> s + getColsCuMat mat) 0 xs
  harmonics <-
    createHarmonicMatriesGPU
      ptxs
      numBatchR2Freqs
      numPoints
      numFreqs
      period
      delta
      deltaFreq
  coefs <- blockMatrixMultiply2 True deviceIDs harmonics xs
  return .
    fromUnboxed (Z :. cols :. numFreqs :. numFreqs) .
    VU.map (* simpsonNorm) . VS.convert . getHostCuMat $
    coefs

computeFourierSeriesR2 ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => [Int]
  -> Int
  -> Int
  -> e
  -> [[CuMat (Complex e)]]
  -> [CuMat (Complex e)]
  -> IO (R.Array U DIM3 (Complex e))
computeFourierSeriesR2 deviceIDs numFreqs numPoints period inverseHarmonics xs = do
  let rows = L.foldl' (\s mat -> s + getRowsCuMat mat) 0 xs
  series <- blockMatrixMultiply1 False deviceIDs xs inverseHarmonics
  return .
    fromUnboxed (Z :. rows :. numPoints :. numPoints) .
    VS.convert . getHostCuMat $
    series


{-# INLINE applyGaussian #-}
applyGaussian ::
     (R.Source s (Complex e), Unbox e, RealFloat e)
  => e
  -> e
  -> R.Array s DIM3 (Complex e)
  -> R.Array D DIM3 (Complex e)
applyGaussian deltaFreq std arr =
  let (Z :. _ :. numFreq :. _) = extent arr
      freq = div numFreq 2
      gaussianCoefficients =
        computeUnboxedS . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
          gaussian2D
            (deltaFreq * (fromIntegral $ xFreq - freq))
            (deltaFreq * (fromIntegral $ yFreq - freq))
            std :+
          0
  in R.traverse2 arr gaussianCoefficients const $ \fC fG idx@(Z :. _ :. i :. j) ->
       fC idx * fG (Z :. i :. j)

-- Stream
{-# INLINE indexSource #-}
indexSource :: Int -> Int -> ConduitT () [(Int, Int)] (ResourceT IO) ()
indexSource numIndex numBatch =
  let xs =
        divideListN numBatch $
        [ (xFreq, yFreq)
        | xFreq <- getListFromNumber numIndex
        , yFreq <- getListFromNumber numIndex
        ]
  in CL.sourceList xs

{-# INLINE fourierCoefficientsConduit #-}
fourierCoefficientsConduit ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , Show e
     )
  => [Int]
  -> [PTX]
  -> Int
  -> e
  -> e
  -> e
  -> [CuMat (Complex e)]
  -> ConduitT [(Int, Int)] (CuMat (Complex e)) (ResourceT IO) ()
fourierCoefficientsConduit deviceIDs ptxs numPoints period delta deltaFreq xs =
  awaitForever $ \idx' -> do
    let idxs = divideListN (L.length deviceIDs) idx'
        harmonics =
          parMap
            rdeepseq
            (uncurry (harmonicMatPTX numPoints period delta deltaFreq)) $
          L.zip ptxs idxs
    liftIO $ printCurrentTime ""
    coefs <- liftIO $ blockMatrixMultiply3 False deviceIDs harmonics xs
    yield $!! coefs

{-# INLINE fourierCoefficientsSink #-}
fourierCoefficientsSink ::
     (Storable e, Unbox e, RealFloat e, Show e)
  => Int
  -> e
  -> ConduitT (CuMat (Complex e)) Void (ResourceT IO) ((R.Array U DIM3 (Complex e)))
fourierCoefficientsSink numR2Freqs delta = do
  xs <- CL.consume
  let cols = getColsCuMat . L.head $ xs
      simpsonNorm = (delta / 3) ^ 2 :+ 0
  return .
    fromUnboxed (Z :. cols :. numR2Freqs :. numR2Freqs) .
    VU.map (* simpsonNorm) .
    VS.convert . getHostCuMat . transposeCuMat . concatCuMat $
    xs

computeFourierCoefficientsR2Stream ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> Int
  -> [CuMat (Complex e)]
  -> IO (R.Array U DIM3 (Complex e))
computeFourierCoefficientsR2Stream deviceIDs ptxs numFreqs numPoints period delta deltaFreq numBatch xs =
  runConduitRes $
  indexSource numFreqs numBatch .|
  fourierCoefficientsConduit deviceIDs ptxs numPoints period delta deltaFreq xs .|
  fourierCoefficientsSink numFreqs delta

{-# INLINE fourierSeriesConduit #-}
fourierSeriesConduit ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , Show e
     )
  => [Int]
  -> [PTX]
  -> Int
  -> e
  -> e
  -> [CuMat (Complex e)]
  -> ConduitT [(Int, Int)] (CuMat (Complex e)) (ResourceT IO) ()
fourierSeriesConduit deviceIDs ptxs numFreqs period delta xs =
  awaitForever $ \idx' -> do
    let idxs = divideListN (L.length deviceIDs) idx'
        harmonics =
          parMap
            rdeepseq
            (uncurry (inverseHarmonicMatPTX1 numFreqs period delta)) $
          L.zip ptxs idxs
    liftIO $ printCurrentTime "fourierSeriesConduit"
    coefs <- liftIO $ blockMatrixMultiply3 False deviceIDs harmonics xs
    yield $!! coefs

{-# INLINE fourierSeriesSink #-}
fourierSeriesSink ::
     (Storable e, Unbox e, RealFloat e, Show e)
  => Int
  -> ConduitT (CuMat (Complex e)) Void (ResourceT IO) ((R.Array U DIM3 (Complex e)))
fourierSeriesSink numPoints = do
  xs <- CL.consume
  let cols = getColsCuMat . L.head $ xs
  return .
    fromUnboxed (Z :. cols :. numPoints :. numPoints) .
    VS.convert . getHostCuMat . transposeCuMat . concatCuMat $
    xs

computeFourierSeriesR2Stream ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> [CuMat (Complex e)]
  -> IO (R.Array U DIM3 (Complex e))
computeFourierSeriesR2Stream deviceIDs ptxs numFreqs numPoints period delta numBatch xs = 
  runConduitRes $
    indexSource numPoints numBatch .|
    fourierSeriesConduit deviceIDs ptxs numFreqs period delta xs .|
    fourierSeriesSink numPoints
    

{-# INLINE fourierSeriesConduitAcc #-}
fourierSeriesConduitAcc ::
     ( A.Floating e
     , A.Elt (Complex e)
     , A.FromIntegral Int e
     , Numeric (Complex e)
     , Unbox e
     )
  => [PTX]
  -> Int
  -> e
  -> e
  -> A.Acc (Matrix (Complex e))
  -> ConduitT [(Int, Int)] (VU.Vector (Complex e)) (ResourceT IO) ()
fourierSeriesConduitAcc ptxs numFreqs period delta x =
  awaitForever $ \idx' -> do
    let idxs = divideListN (L.length ptxs) idx'
        harmonics =
          L.map
            (\r2Positions ->
               inverseHarmonicAcc2 numFreqs period delta .
               A.fromList (A.Z A.:. (L.length r2Positions)) $
               r2Positions) $
          idxs
    liftIO $ printCurrentTime "fourierSeriesConduit"
    yield $!! blockMatrixMultiply ptxs harmonics x

{-# INLINE fourierSeriesSinkAcc #-}
fourierSeriesSinkAcc ::
     (Unbox e)
  => Int
  -> Int
  -> ConduitT (VU.Vector (Complex e)) Void (ResourceT IO) ((R.Array U DIM3 (Complex e)))
fourierSeriesSinkAcc numPoints cols = do
  xs <- CL.consume
  return .
    computeS .
    R.backpermute
      (Z :. cols :. numPoints :. numPoints)
      (\(Z :. a :. b :. c) -> (Z :. b :. c :. a)) .
    fromUnboxed (Z :. numPoints :. numPoints :. cols) . VU.concat $
    xs

computeFourierSeriesR2StreamAcc ::
     ( A.Floating e
     , A.Elt (Complex e)
     , A.FromIntegral Int e
     , Numeric (Complex e)
     , Unbox e
     )
  => [PTX]
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Acc (Matrix (Complex e))
  -> IO (R.Array U DIM3 (Complex e))
computeFourierSeriesR2StreamAcc ptxs numFreqs numPoints cols period delta numBatch x =
  runConduitRes $
  indexSource numPoints numBatch .|
  fourierSeriesConduitAcc ptxs numFreqs period delta (A.compute x) .|
  fourierSeriesSinkAcc numPoints cols

{-# INLINE fourierSeriesConduitAcc' #-}
fourierSeriesConduitAcc' ::
     ( A.Floating e
     , A.Elt (Complex e)
     , A.FromIntegral Int e
     , Numeric (Complex e)
     , Unbox e
     , Prelude.Num e
     )
  => [PTX]
  -> Int
  -> e
  -> e
  -> A.Acc (Matrix (Complex e))
  -> ConduitT [(Int, Int)] (VU.Vector (Complex e)) (ResourceT IO) ()
fourierSeriesConduitAcc' ptxs numFreqs period delta x =
  awaitForever $ \idx' -> do
    let idxs = divideListN (L.length ptxs) idx'
        harmonics =
          L.map
            (\r2Positions ->
               inverseHarmonicAcc2' numFreqs period delta .
               A.fromList (A.Z A.:. L.length r2Positions) $
               r2Positions) $
          idxs
    liftIO $ printCurrentTime "fourierSeriesConduit"
    yield $ blockMatrixMultiply ptxs harmonics x

computeFourierSeriesR2StreamAcc' ::
     ( A.Floating e
     , A.Elt (Complex e)
     , A.FromIntegral Int e
     , Numeric (Complex e)
     , Unbox e
     , Prelude.Num e
     )
  => [PTX]
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Acc (Matrix (Complex e))
  -> IO (R.Array U DIM3 (Complex e))
computeFourierSeriesR2StreamAcc' ptxs numFreqs numPoints cols period delta numBatch x =
  runConduitRes $
  indexSource numPoints numBatch .|
  fourierSeriesConduitAcc' ptxs numFreqs period delta (A.compute x) .|
  fourierSeriesSinkAcc numPoints cols
