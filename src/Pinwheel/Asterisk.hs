{-# LANGUAGE Strict #-}
module Pinwheel.Asterisk where

import Data.Complex

{-# INLINE asteriskGaussianFunc #-}
asteriskGaussianFunc ::
     (RealFloat e) => e -> e -> Int -> e -> e -> e -> Complex e
asteriskGaussianFunc phi rho n a b period =
  (0 :+ (-1)) ^ n *
  exp
    ((-1) *
     (((fromIntegral n ^ 2) / (4 * b) + (pi * rho) ^ 2 / (a * period ^ 2)) :+
      (fromIntegral n * phi)))

