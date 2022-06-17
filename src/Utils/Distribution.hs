{-# LANGUAGE Strict #-}
module Utils.Distribution where

import Data.Complex 

{-# INLINE gaussian1D #-}
gaussian1D :: (Fractional a, Floating a) => a -> a -> a
gaussian1D x std = exp (x ^ 2 / ((-2) * std ^ 2)) / std * sqrt (2 * pi)

{-# INLINE gaussian1DFreq #-}
gaussian1DFreq :: (Fractional a, Floating a) => a -> a -> a
gaussian1DFreq x std = exp ((std * x) ^ 2 / (-2)) -- * std / sqrt (2 * pi)

{-# INLINE gaussian1DFourierCoefficients #-}
gaussian1DFourierCoefficients :: (Fractional a, Floating a) => a -> a -> a 
gaussian1DFourierCoefficients x std =
  std * sqrt (2 * pi) * exp ((std * x) ^ 2 * (-0.5)) 

{-# INLINE gaussian1DLaplaceCoefficients #-}
gaussian1DLaplaceCoefficients :: (RealFloat a) => a -> a -> Complex a 
gaussian1DLaplaceCoefficients x std =
  -- (std * sqrt (2 * pi) :+ 0) *
  exp ((0.5 * std ^ 2 :+ 0) * ((-1) :+ x) ^ 2)

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
