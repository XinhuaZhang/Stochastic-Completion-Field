{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
module FokkerPlanck.Interpolation where

import           Data.Array.Repa as R
import           Data.Complex
import           Data.List       as L
import           Filter.Utils
import           Types

{-# INLINE radialCubicInterpolation #-}
radialCubicInterpolation ::
     (R.Source r1 Double, R.Source r2 (Complex Double), Shape sh)
  => R.Array r1 (sh :. Int) Double
  -> Double
  -> R.Array r2 (sh :. Int :. Int) (Complex Double)
  -> R.Array D (sh :. Int :. Int) (Complex Double)
radialCubicInterpolation radialArr scaleFactor inputArr =
  let (_ :. n) = extent radialArr
      (_ :. nx :. ny) = extent inputArr
   in if scaleFactor == 0
        then error "radialCubicInterpolation: scaleFactor = 0"
        else R.traverse2 inputArr radialArr const $ \fInput fRadial idx@(sh :. i :. j) ->
               let r =
                     (sqrt . fromIntegral $
                      (i - center nx) ^ 2 + (j - center ny) ^ 2) /
                     scaleFactor
                in if r < 0 || r > (fromIntegral $ n - 1) -- out of boundary. When r = 0, log r makes no sense, and therefore the information at r = 0 should not be used to interpoalte 0 < r <= 1.
                     then 0
                     else if r <= 1 -- linear interpolation
                            then fInput idx *
                                 ((fRadial (sh :. 1) +
                                   (r - 1) *
                                   (fRadial (sh :. 2) - fRadial (sh :. 1))) :+
                                  0)
                            else if r >= (fromIntegral $ n - 2) -- linear interpolation
                                   then fInput idx *
                                        ((fRadial (sh :. n - 2) +
                                          (r - (fromIntegral $ n - 2)) *
                                          (fRadial (sh :. n - 1) -
                                           fRadial (sh :. n - 2))) :+
                                         0)
                                   else let x0 = floor r -- cubic interpolation
                                            p0 = fRadial (sh :. x0 - 1)
                                            p1 = fRadial (sh :. x0)
                                            p2 = fRadial (sh :. x0 + 1)
                                            p3 = fRadial (sh :. x0 + 2)
                                            x = r - fromIntegral x0
                                         in fInput idx *
                                            (((-0.5 * p0 + 1.5 * p1 - 1.5 * p2 +
                                               0.5 * p3) *
                                              x ^ 3 +
                                              (p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3) *
                                              x ^ 2 +
                                              (-0.5 * p0 + 0.5 * p2) * x +
                                              p1) :+
                                             0)

{-# INLINE cubicInterpolation #-}
cubicInterpolation ::
     (R.Source s Double)
  => R.Array s DIM5 Double
  -> Double
  -> R.Array D DIM4 Double
cubicInterpolation arr r =
  let (Z :. _ :. _ :. _ :. _ :. n) = extent arr
   in R.traverse arr (\(sh :. _) -> sh) $ \f sh@(Z :. tf :. sf :. t0f :. s0f) ->
        if r < 0 || r > (fromIntegral $ n - 1) -- out of boundary. When r = 0, log r makes no sense, and therefore the information at r = 0 should not be used to interpoalte 0 < r <= 1.
          then 0
          else if r <= 1 -- linear interpolation
                 then (f (sh :. 1) + (r - 1) * (f (sh :. 2) - f (sh :. 1)))
                 else if r >= (fromIntegral $ n - 2) -- linear interpolation
                        then (f (sh :. n - 2) +
                              (r - (fromIntegral $ n - 2)) *
                              (f (sh :. n - 1) - f (sh :. n - 2)))
                        else let x0 = floor r -- cubic interpolation
                                 p0 = f (sh :. x0 - 1)
                                 p1 = f (sh :. x0)
                                 p2 = f (sh :. x0 + 1)
                                 p3 = f (sh :. x0 + 2)
                                 x = r - fromIntegral x0
                              in ((-0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3) *
                                  x ^ 3 +
                                  (p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3) * x ^ 2 +
                                  (-0.5 * p0 + 0.5 * p2) * x +
                                  p1)
