{-# LANGUAGE Strict #-}
module Utils.Distribution where

{-# INLINE gaussian1D #-}
gaussian1D :: (Fractional a, Floating a) => a -> a -> a
gaussian1D x std = (exp (x ^ 2 / ((-2) * std ^ 2))) / (std * sqrt (2 * pi))

{-# INLINE gaussian2D #-}
gaussian2D :: (Fractional a, Floating a) => a -> a -> a -> a
gaussian2D x y std =
  (exp ((x ^ 2 + y ^ 2) / ((-2) * std ^ 2))) / (2 * pi * std ^ 2)
