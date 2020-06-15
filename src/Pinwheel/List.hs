module Pinwheel.List where

import           Data.List as L

{-# INLINE pinwheelFreqsBound #-}
pinwheelFreqsBound :: Int -> Int -> (Int,Int)
pinwheelFreqsBound m n = (-(m + n), m + n)

{-# INLINE pinwheelFreqs #-}
pinwheelFreqs :: Int -> Int -> [Int]
pinwheelFreqs m n = [-(m + n) .. m + n]
