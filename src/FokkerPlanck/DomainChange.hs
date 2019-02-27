{-# LANGUAGE FlexibleContexts #-}
module FokkerPlanck.DomainChange where

import           Data.Array.Repa       as R
import           Data.Complex
import           Data.List             as L
import           FokkerPlanck.Types
import           Numeric.LinearAlgebra as NL
import           Types
import           Utils.Array

{-# INLINE r2z1t0Tor2s1t0 #-}
r2z1t0Tor2s1t0 :: Int -> [Double] -> R2Z1T0Array -> R2S1T0Array
r2z1t0Tor2s1t0 numOrientations freqs (RepaArray arr) =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numFreqs :. numT0Freqs :. xLen :. yLen) = extent arr
      mat1 =
        (numOrientations >< numFreqs) . R.toList $
        R.traverse
          (fromListUnboxed (Z :. numFreqs) freqs)
          (const (Z :. numOrientations :. numFreqs))
          (\f (Z :. i :. j) ->
             exp (0 :+ (-1) * (deltaTheta * fromIntegral i) * f (Z :. j)))
      mat2 = (numFreqs >< (xLen * yLen * numT0Freqs)) . R.toList $ arr
   in fromListUnboxed (Z :. numOrientations :. numT0Freqs :. xLen :. yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2
      
{-# INLINE r2z1Tor2s1 #-}
r2z1Tor2s1 ::
     (Source s (Complex Double))
  => Int
  -> [Double]
  -> R.Array s DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
r2z1Tor2s1 numOrientations freqs arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numFreqs :. xLen :. yLen) = extent arr
      mat1 =
        (numOrientations >< numFreqs) . R.toList $
        R.traverse
          (fromListUnboxed (Z :. numFreqs) freqs)
          (const (Z :. numOrientations :. numFreqs))
          (\f (Z :. i :. j) ->
             exp (0 :+ (-1) * (deltaTheta * fromIntegral i) * f (Z :. j)))
      mat2 = (numFreqs >< (xLen * yLen)) . R.toList $ arr
   in fromListUnboxed (Z :. numOrientations :. xLen :. yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2
