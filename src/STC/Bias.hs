module STC.Bias where

import           Array.UnboxedArray  as AU
import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           DFT.Plan
import           Types

computeBiasR2T0 :: Int -> Int -> [Double] -> [R2S1RPPoint] -> R2T0Array
computeBiasR2T0 xLen yLen theta0Freqs xs =
  let numTheta0Freqs = L.length theta0Freqs
      xShift = div xLen 2
      (xMin, xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      yShift = div yLen 2
      (yMin, yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      vec =
        toUnboxedVector .
        AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
        L.map (\(R2S1RPPoint (x, y, _, _)) -> ((x, y), 1)) $
        xs
   in computeS .
      R.traverse
        (fromUnboxed (Z :. xLen :. yLen) vec)
        (const (Z :. numTheta0Freqs :. xLen :. yLen)) $ \f idx@(Z :. _ :. i :. j) ->
        f (Z :. i :. j)
          
computeBiasR2T0S0 ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> [R2S1RPPoint]
  -> R2T0S0Array
computeBiasR2T0S0 xLen yLen theta0Freqs scale0Freqs xs =
  let numTheta0Freqs = L.length theta0Freqs
      numScale0Freqs = L.length scale0Freqs
      xShift = div xLen 2
      (xMin, xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      yShift = div yLen 2
      (yMin, yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      vec =
        VU.concat .
        L.map
          (\(t0f, s0f) ->
             toUnboxedVector .
             AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
             L.map (\(R2S1RPPoint (x, y, _, _)) -> ((x, y), 1)) $
             xs) $
        [(t0f, s0f) | t0f <- theta0Freqs, s0f <- scale0Freqs]
   in fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) vec
   
   
computeBiasR2T0S0Gaussian ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> [R2S1RPPoint]
  -> R2T0S0Array
computeBiasR2T0S0Gaussian xLen yLen theta0Freqs scale0Freqs theta sigma xs =
  let numTheta0Freqs = L.length theta0Freqs
      numScale0Freqs = L.length scale0Freqs
      xShift = div xLen 2
      (xMin, xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      yShift = div yLen 2
      (yMin, yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      vec =
        VU.concat .
        L.map
          (\(t0f, s0f) ->
             toUnboxedVector .
             AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
             L.map
               (\(R2S1RPPoint (x, y, t, _)) ->
                  ( (x, y)
                  , if s0f == 0
                      then (exp (-2 * pi * t0f ^ 2 * sigma ^ 2) :+ 0) *
                           exp (0 :+ t0f * (theta + t) / 180 * pi)
                      else 0)) $
             xs) $
        [(t0f, s0f) | t0f <- theta0Freqs, s0f <- scale0Freqs]
   in fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) vec
