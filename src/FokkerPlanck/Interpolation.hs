{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
module FokkerPlanck.Interpolation where

import           Data.Array.Repa as R
import           Data.Complex
import           Data.List       as L
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
   in R.traverse2 inputArr radialArr const $ \fInput fRadial idx@(sh :. i :. j) ->
        let r =
              (sqrt . fromIntegral $ (i - div nx 2) ^ 2 + (j - div ny 2) ^ 2) /
              scaleFactor
         in if scaleFactor == 0
              then error "radialCubicInterpolation: scaleFactor = 0"
              else if r < 0 || r > (fromIntegral $ n - 1) -- out of boundary
                     then 0
                     else if r <= 1 -- linear interpolation
                            then fInput idx *
                                 ((fRadial (sh :. 0) +
                                   r * (fRadial (sh :. 1) - fRadial (sh :. 0))) :+
                                  0)
                            else if r >= (fromIntegral $ n - 2) -- linear interpolation
                                   then fInput idx *
                                        ((fRadial (sh :. (n - 2)) +
                                          (r - (fromIntegral $ n - 2)) *
                                          (fRadial (sh :. n - 1) -
                                           fRadial (sh :. n - 2))) :+
                                         0)
                                   else let x0 = floor r -- cubic interpolation
                                            p0 = fRadial (sh :. (x0 - 1))
                                            p1 = fRadial (sh :. x0)
                                            p2 = fRadial (sh :. (x0 + 1))
                                            p3 = fRadial (sh :. (x0 + 2))
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
