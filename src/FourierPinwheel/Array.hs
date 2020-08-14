{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
module FourierPinwheel.Array where

import           Data.Array.Accelerate as A
import           Data.Vector.Generic   as VG
import           Utils.Parallel

data FPArray vector = FPArray
  { getFPArrayNumXFreq     :: Int
  , getFPArrayNumYFreq     :: Int
  , getFPArrayNumRFreq     :: Int
  , getFPArrayNumThetaFreq :: Int
  , getFPArrayNumRhoFreq   :: Int
  , getFPArrayNumPhiFreq   :: Int
  , getFPArray             :: [vector]
  }

{-# INLINE toMatrixAcc #-}
toMatrixAcc ::
     (VG.Vector vector e, Elt e) => FPArray (vector e) -> Acc (A.Array DIM2 e)
toMatrixAcc (FPArray numXFreq numYFreq numRFreq numThetaFreq _ _ arr) =
  A.use .
  A.fromList (Z :. (numRFreq * numThetaFreq) :. (numXFreq * numYFreq)) .
  VG.toList . VG.concat $
  arr

{-# INLINE parMapFPArray #-}
parMapFPArray ::
     (NFData vector) => (vector -> vector) -> FPArray vector -> FPArray vector
parMapFPArray f (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs) =
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
  parMap rdeepseq f $
  vecs
  
{-# INLINE parZipWithFPArray #-}
parZipWithFPArray ::
     (NFData vector)
  => (vector -> vector -> vector)
  -> FPArray vector
  -> FPArray vector
  -> FPArray vector
parZipWithFPArray f (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs1) arr =
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
  parZipWith rdeepseq f vecs1 . getFPArray $
  arr
