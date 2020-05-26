{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Strict              #-}
module Pinwheel.FourierSeries2D where

import           Data.Array.Repa          as R
import           Data.Complex
import           Data.List                as L
import           Data.Vector.Storable     as VS
import           Data.Vector.Unboxed      as VU
import           Foreign.CUDA.Driver      as CUDA
import           Pinwheel.Base
import           Pinwheel.BlockCudaMatrix
import           Utils.BLAS
import           Utils.Distribution
import           Utils.List
import           Utils.Parallel
import           Utils.SimpsonRule
import Image.IO

{-# INLINE harmonicMat #-}
harmonicMat ::
     (Floating e, Unbox e, Storable e)
  => Int
  -> e
  -> e
  -> [(Int, Int)]
  -> CuMat (Complex e)
harmonicMat numPoints period delta r2Freqs =
  let centerPoint = div numPoints 2
      constant = (-2) * pi * delta / period
      rows = L.length r2Freqs
  in CuMat rows (numPoints ^ 2) .
     CuVecHost .
     VU.convert .
     toUnboxed .
     computeUnboxedS .
     R.traverse
       (fromListUnboxed (Z :. rows) r2Freqs)
       (\_ -> (Z :. rows :. numPoints :. numPoints)) $ \f (Z :. r :. x' :. y') ->
       let x = x' - centerPoint
           y = y' - centerPoint
       in cis $
          constant * fromIntegral (fst (f (Z :. r)) * x + snd (f (Z :. r)) * y)

-- create column-major inverse harmonics
{-# INLINE inverseHarmonicMat #-}
inverseHarmonicMat ::
     (Floating e, Unbox e, Storable e)
  => Int
  -> e
  -> e
  -> [(Int, Int)]
  -> CuMat (Complex e)
inverseHarmonicMat numFreqs period delta r2Positions =
  let centerFreq = div numFreqs 2
      constant = 2 * pi * delta / period
      cols = L.length r2Positions
  in CuMat (numFreqs ^ 2) cols .
     CuVecHost .
     VU.convert .
     toUnboxed .
     computeS .
     R.traverse
       (fromListUnboxed (Z :. cols) r2Positions)
       (\_ -> (Z :. numFreqs :. numFreqs :. cols)) $ \f (Z :. xFreq' :. yFreq' :. c) ->
       let xFreq = xFreq' - centerFreq
           yFreq = yFreq' - centerFreq
       in cis $
          constant *
          fromIntegral (fst (f (Z :. c)) * xFreq + snd (f (Z :. c)) * yFreq)

{-# INLINE blockPinwheel #-}
blockPinwheel ::
     forall e. (Storable e, Unbox e, RealFloat e, CUBLAS (Complex e))
  => Int
  -> e
  -> e
  -> [(Int, Int)]
  -> CuMat (Complex e)
blockPinwheel numPoints delta sigma freqs =
  let simpsonWeights =
        computeUnboxedS $
        (computeWeightArrFromListOfShape [numPoints, numPoints] :: R.Array D DIM2 (Complex e))
      numFreqs = L.length freqs
      center = div numPoints 2
      pinwheelArray =
        R.traverse
          (fromListUnboxed (Z :. numFreqs) freqs)
          (\_ -> (Z :. numPoints :. numPoints :. numFreqs)) $ \f (Z :. x :. y :. freq) ->
          fourierMellin
            sigma
            (snd . f $ (Z :. freq))
            (fst . f $ (Z :. freq))
            ( delta * fromIntegral (x - center)
            , delta * fromIntegral (y - center))
      weightedPinwheelArray =
        R.traverse2 pinwheelArray simpsonWeights const $ \fP fS idx@(Z :. x :. y :. _) ->
          fP idx * fS (Z :. x :. y)
  in CuMat (numPoints ^ 2) numFreqs .
     CuVecHost . VU.convert . toUnboxed . computeS $
     weightedPinwheelArray

{-# INLINE applyGaussian #-}
applyGaussian ::
     (R.Source s (Complex e), Unbox e, RealFloat e)
  => e
  -> R.Array s DIM4 (Complex e)
  -> IO (R.Array U DIM4 (Complex e))
applyGaussian std arr =
  let (Z :. _ :. _ :. numFreq :. _) = extent arr
      freq = div numFreq 2
      gaussianCoefficients =
        computeUnboxedS . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
          gaussian2D
            (fromIntegral $ xFreq - freq)
            (fromIntegral $ yFreq - freq)
            std :+
          0
  in computeUnboxedP . R.traverse2 arr gaussianCoefficients const $ \fC fG idx@(Z :. _ :. _ :. i :. j) ->
       fC idx * fG (Z :. i :. j)

-- The input xs is column-major
computeFourierCoefficients ::
     forall e. (Storable e, CUBLAS (Complex e), Unbox e, RealFloat e)
  => [Int]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Int
  -> Int
  -> Int
  -> IO (R.Array U DIM4 (Complex e))
computeFourierCoefficients deviceIDs numFreqs numPoints period delta batchSizeR2Freqs batchSizePolarFreqs numAngularFreqs numRadialFreqs = do
  let simpsonNorm = (delta / 3) ^ 2 :+ 0
      harmonics =
        parMap rdeepseq (harmonicMat numPoints period delta) .
        divideList batchSizeR2Freqs $
        [ (xFreq, yFreq)
        | xFreq <- getListFromNumber numFreqs
        , yFreq <- getListFromNumber numFreqs
        ]
      pinwheels =
        parMap
          rdeepseq
          (L.map (blockPinwheel numPoints delta 0.5) .
           divideList batchSizePolarFreqs) .
        divideList
          (ceiling $
           (fromIntegral $ numRadialFreqs * numAngularFreqs) /
           (fromIntegral $ L.length deviceIDs)) $
        [ (radialFreq, angularFreq)
        | radialFreq <- getListFromNumber numRadialFreqs
        , angularFreq <- getListFromNumber numAngularFreqs
        ]
  coefs <- blockMatrixMultiply True deviceIDs harmonics pinwheels
  return .
    fromUnboxed (Z :. numRadialFreqs :. numAngularFreqs :. numFreqs :. numFreqs) .
    VU.map (* simpsonNorm) . VS.convert . getHostCuMat $
    coefs

computeFourierSeries ::
     forall e. (Storable e, CUBLAS (Complex e), Unbox e, RealFloat e)
  => [Int]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Int
  -> Int
  -> Int
  -> [VS.Vector (Complex e)]
  -> IO (R.Array U DIM4 (Complex e))
computeFourierSeries deviceIDs numFreqs numPoints period delta batchSizePolarFreqs batchSizeR2 numAngularFreqs numRadialFreqs coefficients = do
  let centerFreq = div numFreqs 2
      matAs =
        parMap rdeepseq createCuMat . divideList batchSizePolarFreqs $
        coefficients
      invHarmonics =
        parMap
          rdeepseq
          (L.map (inverseHarmonicMat numFreqs period delta) .
           divideList batchSizeR2) .
        divideList
          (ceiling $
           (fromIntegral $ numPoints ^ 2) / (fromIntegral $ L.length deviceIDs)) $
        [ (x, y)
        | x <- getListFromNumber numPoints
        , y <- getListFromNumber numPoints
        ]
  recon <- blockMatrixMultiply False deviceIDs matAs invHarmonics
  return .
    fromUnboxed
      (Z :. numRadialFreqs :. numAngularFreqs :. numPoints :. numPoints) .
    VS.convert . getHostCuMat $
    recon 


-- Utilities
{-# INLINE getNumFreqs #-}
getNumFreqs :: Int -> Int -> Int
getNumFreqs m n = 2 * (m + n) + 1

{-# INLINE getPinwheelFreqs #-}
getPinwheelFreqs :: Int -> Int -> [Int]
getPinwheelFreqs m n =
  let x = m + n
  in [-x .. x]

{-# INLINE getListFromNumber #-}
getListFromNumber :: Int -> [Int]
getListFromNumber n =
  let m = div n 2
  in if odd n
       then [-m .. m]
       else [-m .. (m - 1)]

printCuMat ::
     (RealFloat e, Storable e, Unbox e)
  => FilePath
  -> Bool
  -> Int
  -> CuMat (Complex e)
  -> IO ()
printCuMat filePath isRow m (CuMat rows cols (CuVecHost vec)) =
  let n =
        round . sqrt . fromIntegral $
        if isRow
          then cols
          else rows
      arr = fromUnboxed (Z :. rows :. cols) . VS.convert $ vec
  in plotImageRepaComplex filePath .
     ImageRepa 8 . computeS . R.reshape (Z :. (1 :: Int) :. n :. n) $
     if isRow
       then R.slice arr (Z :. m :. All)
       else R.slice arr (Z :. All :. m)
