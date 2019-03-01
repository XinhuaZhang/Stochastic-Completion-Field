module Array.UnboxedArray where

import           Data.Ix
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           GHC.Arr             as Arr

data UnboxedArray i e =
  UnboxedArray !i
               !i
               (Vector e)
instance (Ix a1, Ix a2, Ix a3, Ix a4, Ix a5, Ix a6) =>
         Ix (a1, a2, a3, a4, a5, a6) where
  range ((l1, l2, l3, l4, l5, l6), (u1, u2, u3, u4, u5, u6)) =
    [ (i1, i2, i3, i4, i5, i6)
    | i1 <- range (l1, u1)
    , i2 <- range (l2, u2)
    , i3 <- range (l3, u3)
    , i4 <- range (l4, u4)
    , i5 <- range (l5, u5)
    , i6 <- range (l6, u6)
    ]
  unsafeIndex ((l1, l2, l3, l4, l5, l6), (u1, u2, u3, u4, u5, u6)) (i1, i2, i3, i4, i5, i6) =
    Arr.unsafeIndex (l6, u6) i6 +
    unsafeRangeSize (l6, u6) *
    (Arr.unsafeIndex (l5, u5) i5 +
     unsafeRangeSize (l5, u5) *
     (Arr.unsafeIndex (l4, u4) i4 +
      unsafeRangeSize (l4, u4) *
      (Arr.unsafeIndex (l3, u3) i3 +
       unsafeRangeSize (l3, u3) *
       (Arr.unsafeIndex (l2, u2) i2 +
        unsafeRangeSize (l2, u2) * (Arr.unsafeIndex (l1, u1) i1)))))
  inRange ((l1, l2, l3, l4, l5, l6), (u1, u2, u3, u4, u5, u6)) (i1, i2, i3, i4, i5, i6) =
    inRange (l1, u1) i1 &&
    inRange (l2, u2) i2 &&
    inRange (l3, u3) i3 &&
    inRange (l4, u4) i4 && inRange (l5, u5) i5 && inRange (l6, u6) i6

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
