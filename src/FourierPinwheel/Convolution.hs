{-# LANGUAGE Strict #-}
module FourierPinwheel.Convolution where

import           Data.Array.Repa           as R
import           Data.List                 as L
import           Data.Vector.Storable      as VS
import           Data.Vector.Unboxed       as VU
import           FourierPinwheel.Array
import           FourierPinwheel.Harmonics
import           Utils.BLAS
import           Utils.Parallel

-- Convolution using BLAS
convolve ::
     (Storable e, Num e, BLAS e, Unbox e)
  => VS.Vector e
  -> FPHarmonicsVector VS.Vector e
  -> FPArray (VS.Vector e)
  -> IO (FPArray (VS.Vector e))
convolve coefficients (FPHarmonics harmonics harmonicsOffset harmonicsNorm) (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq inputArr) = do
  let inputVecs =
        parZipWith rdeepseq (VS.zipWith (*)) harmonics $
        if numRFreq == 1
          then L.concat . L.replicate numRhoFreq $ inputArr
          else inputArr
      m = numRFreq * numThetaFreq
      n = numXFreq * numYFreq
      k = numRhoFreq * numPhiFreq
  vec2 <- gemmBLAS m n k coefficients . VS.concat $ inputVecs
  let arr =
        fromUnboxed (Z :. numRFreq :. numThetaFreq :. numXFreq :. numYFreq) .
        VS.convert $
        vec2
      outputArr =
        parZipWith3
          rdeepseq
          (\(i, j) offsetVec norm ->
             VS.map (\x -> x - norm) .
             VS.zipWith (*) offsetVec .
             VS.convert . toUnboxed . computeS . R.slice arr $
             (Z :. i :. j :. All :. All))
          [(i, j) | i <- [0 .. numRFreq - 1], j <- [0 .. numThetaFreq - 1]]
          harmonicsOffset
          harmonicsNorm
  return . FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq $
    outputArr
