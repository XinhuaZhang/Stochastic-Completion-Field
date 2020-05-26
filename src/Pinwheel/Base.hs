{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Strict #-}
module Pinwheel.Base where

import           Data.Complex

{-# INLINE fourierMellin #-}
fourierMellin :: (Eq e, RealFloat e) => e -> Int -> Int -> (e, e) -> Complex e
fourierMellin sigma angularFreq radialFreq (!x, !y) =
  if x == 0 && y == 0
    then 0
    else (x :+ y) ** (fromIntegral angularFreq :+ 0) *
         ((x ^ 2 + y ^ 2) :+ 0) **
         (((-(fromIntegral angularFreq) - 1 + sigma) :+ fromIntegral radialFreq) /
          2)


