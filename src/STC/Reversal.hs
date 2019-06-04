{-# LANGUAGE FlexibleContexts #-}
module STC.Reversal where

import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Types

{-# INLINE computeReversalR2Z2T0S0 #-}
computeReversalR2Z2T0S0 ::
     (R.Source r (Complex Double))
  => [Double]
  -> R.Array r DIM6 (Complex Double)
  -> R.Array D DIM6 (Complex Double)
computeReversalR2Z2T0S0 theta0Freqs sourceArr =
  R.traverse2
    sourceArr
    (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
    const $ \f1 f2 idx@(Z :. _ :. _ :. l :. _ :. i :. j) ->
    f1 idx * exp (0 :+ f2 (Z :. l) * pi)
