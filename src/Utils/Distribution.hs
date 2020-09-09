{-# LANGUAGE Strict #-}
module Utils.Distribution where

{-# INLINE gaussian1D #-}
gaussian1D :: (Fractional a, Floating a) => a -> a -> a
gaussian1D x std = exp (x ^ 2 / ((-2) * std ^ 2)) / std * sqrt (2 * pi)

{-# INLINE gaussian1DFreq #-}
gaussian1DFreq :: (Fractional a, Floating a) => a -> a -> a
gaussian1DFreq x std = exp ((std * x) ^ 2 / (-2)) * std / sqrt (2 * pi)

{-# INLINE gaussian1DFourierCoefficients #-}
gaussian1DFourierCoefficients :: (Fractional a, Floating a) => a -> a -> a -> a
gaussian1DFourierCoefficients x p std =
  exp ((pi * std * x / p) ^ 2 * (-2)) * std * sqrt (2 * pi) / p

{-# INLINE gaussian2D #-}
gaussian2D :: (Fractional a, Floating a) => a -> a -> a -> a
gaussian2D x y std = exp ((x ^ 2 + y ^ 2) / ((-2) * std ^ 2)) / 2 * pi * std ^ 2

{-# INLINE gaussian2DPolar #-}
gaussian2DPolar :: (Fractional a, Floating a) => a -> a -> a
gaussian2DPolar rho std = exp (rho ^ 2 / ((-2) * std ^ 2)) 

{-# INLINE gaussian2DFourierCoefficients #-}
gaussian2DFourierCoefficients ::
     (Fractional a, Floating a) => a -> a -> a -> a -> a
gaussian2DFourierCoefficients x y p std =
  exp ((x ^ 2 + y ^ 2) * (-2) * (pi * std / p) ^ 2) * 2 * pi * std ^ 2 / p ^ 2
