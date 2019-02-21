module Array.UnboxedArray where

import           Data.Ix
import           Data.List           as L
import           Data.Vector.Unboxed as VU

data UnboxedArray i e =
  UnboxedArray !i
               !i
               (Vector e)

{-# INLINE accum #-}
accum ::
     (Ix i, Unbox a, Unbox b)
  => (a -> b -> a)
  -> a
  -> (i, i)
  -> [(i, b)]
  -> UnboxedArray i a
accum f init range'@(l, u) =
  UnboxedArray l u .
  VU.accum f (VU.replicate (rangeSize range') init) .
  L.map (\(idx, x) -> (index range' idx, x))

{-# INLINE accumulate #-}
accumulate ::
     (Ix i, Unbox i, Unbox a, Unbox b)
  => (a -> b -> a)
  -> a
  -> (i, i)
  -> Vector (i, b)
  -> UnboxedArray i a
accumulate f init range'@(l, u) =
  UnboxedArray l u .
  VU.accumulate f (VU.replicate (rangeSize range') init) .
  VU.map (\(idx, x) -> (index range' idx, x))

{-# INLINE toUnboxedVector #-}
toUnboxedVector :: (Ix i, Unbox a) => UnboxedArray i a -> Vector a
toUnboxedVector (UnboxedArray _ _ vec) = vec
