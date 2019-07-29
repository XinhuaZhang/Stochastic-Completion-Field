{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module FokkerPlanck.Pinwheel where

import           Data.Array.Repa            as R
import           Data.Complex
import           Data.List                  as L
import           Filter.Pinwheel
import           FokkerPlanck.Interpolation
import           Types
import           Utils.Coordinates


{-# INLINE computeR2Z1T0ArrayRadial #-}
computeR2Z1T0ArrayRadial ::
     (R.Source r Double)
  => R.Array r DIM3 Double
  -> Int
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> R.Array D DIM4 (Complex Double)
computeR2Z1T0ArrayRadial radialArr xLen yLen scaleFactor thetaFreqs theta0Freqs =
  let pinwheelArr =
        traverse2
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
          (\(Z :. numThetaFreq) (Z :. numTheta0Freq) ->
             (Z :. numThetaFreq :. numTheta0Freq :. xLen :. yLen)) $ \f f0 (Z :. k :. l :. i :. j) ->
          pinwheel
            (f (Z :. k) - f0 (Z :. l))
            0
            (exp 1)
            0
            (i - center xLen)
            (j - center yLen)
   in radialCubicInterpolation radialArr scaleFactor pinwheelArr


{-# INLINE computeR2Z2T0S0ArrayRadial #-}
computeR2Z2T0S0ArrayRadial ::
     (R.Source r Double)
  => R.Array r DIM5 Double
  -> Int
  -> Int
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> R.Array D DIM6 (Complex Double)
computeR2Z2T0S0ArrayRadial radialArr xLen yLen scaleFactor rMax thetaFreqs scaleFreqs theta0Freqs scale0Freqs = do
  let pinwheelArr =
        traverse4
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          (fromListUnboxed (Z :. L.length scaleFreqs) scaleFreqs)
          (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
          (fromListUnboxed (Z :. L.length scale0Freqs) scale0Freqs)
          (\(Z :. numThetaFreq) (Z :. numScaleFreq) (Z :. numTheta0Freq) (Z :. numScale0Freq) ->
             (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :.
              numScale0Freq :.
              xLen :.
              yLen)) $ \ft fs ft0 fs0 (Z :. t :. s :. t0 :. s0 :. i :. j) ->
          pinwheel
            (ft (Z :. t) - ft0 (Z :. t0))
            (fs (Z :. s) + fs0 (Z :. s0))
            rMax
            0
            (i - center xLen)
            (j - center yLen)
   in radialCubicInterpolation radialArr scaleFactor pinwheelArr
   
{-# INLINE computeR2Z2T0S0ArrayRadial' #-}
computeR2Z2T0S0ArrayRadial' ::
     Int
  -> Int
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> R.Array D DIM6 (Complex Double)
computeR2Z2T0S0ArrayRadial' xLen yLen scaleFactor rMax thetaFreqs scaleFreqs theta0Freqs scale0Freqs = do
  let pinwheelArr =
        traverse4
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          (fromListUnboxed (Z :. L.length scaleFreqs) scaleFreqs)
          (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
          (fromListUnboxed (Z :. L.length scale0Freqs) scale0Freqs)
          (\(Z :. numThetaFreq) (Z :. numScaleFreq) (Z :. numTheta0Freq) (Z :. numScale0Freq) ->
             (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :.
              numScale0Freq :.
              xLen :.
              yLen)) $ \ft fs ft0 fs0 (Z :. t :. s :. t0 :. s0 :. i :. j) ->
          pinwheel
            (ft (Z :. t) - ft0 (Z :. t0))
            (fs (Z :. s) - fs0 (Z :. s0))
            rMax
            0
            (i - center xLen)
            (j - center yLen)
   in pinwheelArr

{-# INLINE cutoff #-}
cutoff :: (R.Source r Double) => Int -> R.Array r DIM5 Double -> R.Array D DIM5 Double
cutoff r arr =
  R.traverse arr id $ \f idx@(Z :. _ :. _ :. _ :. _ :. e) ->
    if e > r 
      then 0
      else f idx
