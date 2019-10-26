{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
module FokkerPlanck.Pinwheel
  ( module Filter.Pinwheel
  , module FokkerPlanck.Pinwheel
  ) where

import           Data.Array.Repa               as R
import           Data.Complex
import           Data.List                     as L
import           Data.Vector.Unboxed           as VU
import           Filter.Pinwheel
import           FokkerPlanck.Interpolation
import           Numeric.LinearAlgebra.Data    as NL
import           Numeric.LinearAlgebra.HMatrix
import           Types
import           Utils.Array
import           Utils.Coordinates
import           Utils.Parallel

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
  => PinwheelType
  -> R.Array r DIM5 Double
  -> Int
  -> Int
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> R.Array D DIM6 (Complex Double)
computeR2Z2T0S0ArrayRadial pinwheelType radialArr xLen yLen scaleFactor rMax thetaFreqs scaleFreqs theta0Freqs scale0Freqs =
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
          pinwheelFunc
            pinwheelType
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
computeR2Z2T0S0ArrayRadial' xLen yLen scaleFactor rMax thetaFreqs scaleFreqs theta0Freqs scale0Freqs =
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

{-# INLINE computeLocalEigenVector #-}
computeLocalEigenVector ::
     (R.Source r Double)
  => ParallelParams
  -> (Double -> Double -> Double -> Double -> Int -> Int -> Complex Double)
  -> R.Array r DIM5 Double
  -> Int
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
computeLocalEigenVector parallelParams pinwheelFunc radialArr xLen yLen rMax thetaFreqs scaleFreqs =
  let (Z :. numThetaFreq :. numScaleFreq :. _ :. _ :. _) = extent radialArr
      len = numThetaFreq * numScaleFreq
   in computeS .
      rotate4D .
      rotate4D .
      fromUnboxed (Z :. xLen :. yLen :. numThetaFreq :. numScaleFreq) .
      VU.concat .
      parMapChunk
        parallelParams
        rdeepseq
        (\(x, y) ->
           let r =
                 sqrt $
                 (fromIntegral $ x - center xLen) ^ 2 +
                 (fromIntegral $ y - center yLen) ^ 2
               interpolatedArray = cubicInterpolation radialArr r
               interpolatedPinwheel =
                 R.traverse3
                   interpolatedArray
                   (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
                   (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
                   (\sh _ _ -> sh) $ \f ft fs idx@(Z :. t :. s :. t0 :. s0) ->
                   (f idx :+ 0) *
                   pinwheelFunc
                     (ft (Z :. t) - ft (Z :. t0))
                     (fs (Z :. s) + fs (Z :. s0))
                     rMax
                     0
                     (x - center xLen)
                     (y - center yLen)
               (eigVal, eigVec) =
                 eig . (len >< len) . R.toList $ interpolatedPinwheel
               (val, vec) =
                 L.head . L.reverse . L.sortOn (magnitude . fst) $
                 L.zip (NL.toList eigVal) (toColumns $ eigVec)
               dominantVec = VU.fromList . L.map (* val) . NL.toList $ vec
            in -- VU.zipWith
               --   (*)
               --   (VU.fromList
               --      [ exp
               --        (0 :+
               --         (tf) *
               --         (angleFunctionRad
               --            (fromIntegral $ x - center xLen)
               --            (fromIntegral $ y - center yLen)))
               --      | tf <- thetaFreqs
               --      , sf <- scaleFreqs
               --      ])
               --   dominantVec
               dominantVec
        ) $
      [(x, y) | x <- [0 .. xLen - 1], y <- [0 .. yLen - 1]]
      
{-# INLINE computeLocalEigenVectorSink #-}
computeLocalEigenVectorSink ::
     (R.Source r Double)
  => ParallelParams
  -> (Double -> Double -> Double -> Double -> Int -> Int -> Complex Double)
  -> R.Array r DIM5 Double
  -> Int
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
computeLocalEigenVectorSink parallelParams pinwheelFunc radialArr xLen yLen rMax thetaFreqs scaleFreqs =
  let (Z :. numThetaFreq :. numScaleFreq :. _ :. _ :. _) = extent radialArr
      len = numThetaFreq * numScaleFreq
   in computeS .
      rotate4D .
      rotate4D .
      fromUnboxed (Z :. xLen :. yLen :. numThetaFreq :. numScaleFreq) .
      VU.concat .
      parMapChunk
        parallelParams
        rdeepseq
        (\(x, y) ->
           let r =
                 sqrt $
                 (fromIntegral $ x - center xLen) ^ 2 +
                 (fromIntegral $ y - center yLen) ^ 2
               interpolatedArray = cubicInterpolation radialArr r
               interpolatedPinwheel =
                 R.traverse3
                   interpolatedArray
                   (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
                   (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
                   (\sh _ _ -> sh) $ \f ft fs idx@(Z :. t :. s :. t0 :. s0) ->
                   (f idx :+ 0) * (exp (0 :+ (ft (Z :. t) - ft (Z :. t0)) * pi)) *
                   pinwheelFunc 
                     (ft (Z :. t) - ft (Z :. t0))
                     (fs (Z :. s) + fs (Z :. s0))
                     rMax
                     0
                     (x - center xLen)
                     (y - center yLen)
               (eigVal, eigVec) =
                 eig . (len >< len) . R.toList $ interpolatedPinwheel
               (val, vec) =
                 L.head . L.reverse . L.sortOn (realPart . fst) $
                 L.zip (NL.toList eigVal) (toColumns $ eigVec)
            in VU.fromList . L.map (* val) . NL.toList $ vec) $
      [(x, y) | x <- [0 .. xLen - 1], y <- [0 .. yLen - 1]]
