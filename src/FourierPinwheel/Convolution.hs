{-# LANGUAGE Strict #-}
module FourierPinwheel.Convolution where

import           Control.Concurrent.Async
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
  => FPData (VS.Vector e)
  -> FPArray (VS.Vector e)
  -> IO (FPArray (VS.Vector e))
convolve (FPData harmonics harmonicsOffset coefs coefHollows) (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq inputArr) = do
  let inputVecs =
        if numRFreq == 1
          then L.concat . L.replicate numRhoFreq $ inputArr
          else inputArr
      m = numThetaFreq
      n = numXFreq * numYFreq
      k = numRhoFreq * numPhiFreq
      inputVec1 =
        VS.concat . parZipWith rdeepseq (VS.zipWith (*)) harmonics $ inputVecs
      inputVec2 = VS.concat inputVecs
  vecs1 <- mapConcurrently (\coef -> gemmBLAS m n k coef inputVec1) coefs
  vecs2 <-
    mapConcurrently
      (\coefHollow -> gemmBLAS m n k coefHollow inputVec2)
      coefHollows
  let outputArr =
        parZipWith3
          rdeepseq
          (\vec1 vec2 offsetVec ->
             VS.zipWith (-) (VS.zipWith (*) offsetVec vec1) vec2)
          vecs1
          vecs2
          harmonicsOffset
  return . FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq $
    outputArr
