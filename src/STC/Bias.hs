{-# LANGUAGE BangPatterns #-}
module STC.Bias where

import           Array.UnboxedArray   as AU
import           Data.List            as L
import           Data.Vector.Storable as VS
import           Data.Vector.Unboxed  as VU
import           STC.Point
import           STC.Utils

{-# INLINE computeBias #-}
computeBias ::
     (Num a, Unbox a, Storable a) => Int -> Int -> [Point] -> VS.Vector a
computeBias !rows !cols =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
  in VS.convert .
     toUnboxedVector .
     AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
     L.map (\(Point x y theta scale) -> ((round x, round y), 1))
