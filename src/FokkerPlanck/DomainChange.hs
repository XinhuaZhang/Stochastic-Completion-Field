{-# LANGUAGE FlexibleContexts #-}
module FokkerPlanck.DomainChange where

import           Data.Array.Repa       as R
import           Data.Complex
import           Data.List             as L
import           FokkerPlanck.Types
import           Numeric.LinearAlgebra as NL
import           Types
import           Utils.Array

{-# INLINE r2z1t0Tor2s1t0 #-}
r2z1t0Tor2s1t0 :: Int -> [Double] -> R2Z1T0Array -> R2S1T0Array
r2z1t0Tor2s1t0 numOrientations freqs arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numFreqs :. numT0Freqs :. xLen :. yLen) = extent arr
      mat1 =
        (numOrientations >< numFreqs) . R.toList $
        R.traverse
          (fromListUnboxed (Z :. numFreqs) freqs)
          (const (Z :. numOrientations :. numFreqs))
          (\f (Z :. i :. j) ->
             exp (0 :+ (-1) * (deltaTheta * fromIntegral i) * f (Z :. j)))
      mat2 = (numFreqs >< (xLen * yLen * numT0Freqs)) . R.toList $ arr
   in fromListUnboxed (Z :. numOrientations :. numT0Freqs :. xLen :. yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2
      
{-# INLINE r2z1Tor2s1 #-}
r2z1Tor2s1 ::
     (Source s (Complex Double))
  => Int
  -> [Double]
  -> R.Array s DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
r2z1Tor2s1 numOrientations freqs arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numFreqs :. xLen :. yLen) = extent arr
      mat1 =
        (numOrientations >< numFreqs) . R.toList $
        R.traverse
          (fromListUnboxed (Z :. numFreqs) freqs)
          (const (Z :. numOrientations :. numFreqs))
          (\f (Z :. i :. j) ->
             exp (0 :+ (-1) * (deltaTheta * fromIntegral i) * f (Z :. j)))
      mat2 = (numFreqs >< (xLen * yLen)) . R.toList $ arr
   in fromListUnboxed (Z :. numOrientations :. xLen :. yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2
      

{-# INLINE r2s1tor2z1 #-}
r2s1tor2z1 ::
     (Source s (Complex Double))
  => [Double]
  -> R.Array s DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
r2s1tor2z1 freqs arr =
  let (Z :. numOrientations :. xLen :. yLen) = extent arr
      deltaTheta = 2 * pi / (fromIntegral numOrientations)
      numFreqs = L.length freqs
      mat1 =
        (numFreqs >< numOrientations) . R.toList $
        R.traverse
          (fromListUnboxed (Z :. numFreqs) freqs)
          (const (Z :. numFreqs :. numOrientations))
          (\f (Z :. i :. j) ->
             exp (0 :+ (deltaTheta * fromIntegral j) * f (Z :. i)))
      mat2 = (numOrientations >< (xLen * yLen)) . R.toList $ arr
   in fromListUnboxed (Z :. numFreqs :. xLen :. yLen) . NL.toList . flatten $
      mat1 NL.<> mat2


{-# INLINE r2z2t0s0Tor2s1rpt0s0 #-}
r2z2t0s0Tor2s1rpt0s0 ::
     Int -> [Double] -> Int -> [Double] -> R2Z2T0S0Array -> R2S1RPT0S0Array
r2z2t0s0Tor2s1rpt0s0 numOrientations thetaFreqs numScales scaleFreqs arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. xLen :. yLen) =
        extent arr
      mat1 =
        ((numOrientations * numScales) >< (numThetaFreq * numScaleFreq)) .
        R.toList $
        R.traverse3
          (fromFunction
             (Z :. numOrientations :. numScales :. numThetaFreq :. numScaleFreq)
             (const 0))
          (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
          (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
          (\a _ _ -> a)
          (\f1 f2 f3 (Z :. i :. j :. k :. l) ->
             exp
               (0 :+
                (-1) *
                ((deltaTheta * fromIntegral i) * f2 (Z :. k) +
                 (2 * pi * fromIntegral j / fromIntegral numScales) *
                 f3 (Z :. l))))
      mat2 =
        ((numThetaFreq * numScaleFreq) ><
         (xLen * yLen * numTheta0Freq * numScale0Freq)) .
        R.toList $
        arr
   in fromListUnboxed
        (Z :. numOrientations :. numScales :. numTheta0Freq :. numScale0Freq :.
         xLen :.
         yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2


{-# INLINE r2z2Tor2s1rp #-}
r2z2Tor2s1rp ::
     (Source s (Complex Double))
  => Int
  -> [Double]
  -> Int
  -> [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
r2z2Tor2s1rp numOrientations thetaFreqs numScales scaleFreqs arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numThetaFreq :. numScaleFreq :. xLen :. yLen) = extent arr
      mat1 =
        ((numOrientations * numScales) >< (numThetaFreq * numScaleFreq)) .
        R.toList $
        R.traverse3
          (fromFunction
             (Z :. numOrientations :. numScales :. numThetaFreq :. numScaleFreq)
             (const 0))
          (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
          (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
          (\a _ _ -> a)
          (\f1 f2 f3 (Z :. i :. j :. k :. l) ->
             exp
               (0 :+
                (1) *
                ((deltaTheta * fromIntegral i) * f2 (Z :. k) +
                 (2 * pi * fromIntegral j / fromIntegral numScales) *
                 f3 (Z :. l))))
      mat2 = ((numThetaFreq * numScaleFreq) >< (xLen * yLen)) . R.toList $ arr
   in fromListUnboxed (Z :. numOrientations :. numScales :. xLen :. yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2
      

{-# INLINE r2z2Tor2s1rpP #-}
r2z2Tor2s1rpP ::
     (Source s (Complex Double))
  => Int
  -> [Double]
  -> Int
  -> [Double]
  -> Double
  -> R.Array s DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
r2z2Tor2s1rpP numOrientations thetaFreqs numScales scaleFreqs maxR arr = do
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      (Z :. numThetaFreq :. numScaleFreq :. xLen :. yLen) = extent arr
      mat2 = ((numThetaFreq * numScaleFreq) >< (xLen * yLen)) . R.toList $ arr
  mat1 <-
    fmap
      (((numOrientations * numScales) >< (numThetaFreq * numScaleFreq)) .
       R.toList) .
    computeUnboxedP $
    R.traverse3
      (fromFunction
         (Z :. numOrientations :. numScales :. numThetaFreq :. numScaleFreq)
         (const 0))
      (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
      (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
      (\a _ _ -> a)
      (\f1 f2 f3 (Z :. i :. j :. k :. l) ->
         exp
           (0 :+
            (1) *
            ((deltaTheta * fromIntegral i) * f2 (Z :. k) +
             (2 * pi * fromIntegral j / fromIntegral numScales) *
             f3 (Z :. l))))
  return .
    fromListUnboxed (Z :. numOrientations :. numScales :. xLen :. yLen) .
    NL.toList . flatten $
    mat1 NL.<> mat2



{-# INLINE r2z2t0s0Tor2s1rps1rp #-}
r2z2t0s0Tor2s1rps1rp ::
     Int -> [Double] -> [Double] -> Int -> [Double] -> [Double] -> Double -> R2Z2T0S0Array -> R.Array U DIM6 (Complex Double)
r2z2t0s0Tor2s1rps1rp numOrientations thetaFreqs theta0Freqs numScales scaleFreqs scale0Freqs maxScale arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      deltaScale =
        if maxScale == 0
          then 0
          else (2 * pi / (log maxScale)) / (fromIntegral numScales)
      (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. xLen :. yLen) =
        extent arr
      mat1 =
        ((numOrientations * numScales * numOrientations * numScales) ><
         (numThetaFreq * numScaleFreq * numTheta0Freq * numScale0Freq)) .
        R.toList $
        R.traverse4
          (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
          (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
          (fromListUnboxed (Z :. numTheta0Freq) theta0Freqs)
          (fromListUnboxed (Z :. numScale0Freq) scale0Freqs)
          (\_ _ _ _ ->
             (Z :. numOrientations :. numScales :. numOrientations :. numScales :.
              numThetaFreq :.
              numScaleFreq :.
              numTheta0Freq :.
              numScale0Freq))
          (\ft fs ft0 fs0 (Z :. o :. s :. o0 :. s0 :. tf :. sf :. tf0 :. sf0) ->
             exp
               (0 :+
                (-1) *
                ((deltaTheta * fromIntegral o) * ft (Z :. tf) +
                 (deltaTheta * fromIntegral o0) * ft0 (Z :. tf0) +
                 (deltaScale * fromIntegral s) * fs (Z :. sf) +
                 (deltaScale * fromIntegral s0) * fs0 (Z :. sf0))))
      mat2 =
        ((numThetaFreq * numScaleFreq * numTheta0Freq * numScale0Freq) ><
         (xLen * yLen)) .
        R.toList $
        arr
   in fromListUnboxed
        (Z :. numOrientations :. numScales :. numOrientations :. numScales :.
         xLen :.
         yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2


{-# INLINE r2s1rps1rpTor2z2t0s0 #-}
r2s1rps1rpTor2z2t0s0 ::
     Int
  -> [Double]
  -> [Double]
  -> Int
  -> [Double]
  -> [Double]
  -> Double
  -> R.Array U DIM6 (Complex Double)
  -> R2Z2T0S0Array
r2s1rps1rpTor2z2t0s0 numOrientations thetaFreqs theta0Freqs numScales scaleFreqs scale0Freqs maxScale arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      deltaScale =
        if maxScale == 0
          then 0
          else (2 * pi / (log maxScale)) / (fromIntegral numScales)
      numThetaFreq = L.length thetaFreqs
      numScaleFreq = L.length scaleFreqs
      numTheta0Freq = L.length theta0Freqs
      numScale0Freq = L.length scale0Freqs
      (Z :. _ :. _ :. _ :. _ :. xLen :. yLen) =
        extent arr
      mat1 =
        ((numThetaFreq * numScaleFreq * numTheta0Freq * numScale0Freq) ><
         (numOrientations * numScales * numOrientations * numScales)) .
        R.toList $
        R.traverse4
          (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
          (fromListUnboxed (Z :. numScaleFreq) scaleFreqs)
          (fromListUnboxed (Z :. numTheta0Freq) theta0Freqs)
          (fromListUnboxed (Z :. numScale0Freq) scale0Freqs)
          (\_ _ _ _ ->
             (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :.
              numScale0Freq :.
              numOrientations :.
              numScales :.
              numOrientations :.
              numScales))
          (\ft fs ft0 fs0 (Z :. tf :. sf :. tf0 :. sf0 :. o :. s :. o0 :. s0) ->
             exp
               (0 :+ (1) *
                ((deltaTheta * fromIntegral o) * ft (Z :. tf) +
                 (deltaTheta * fromIntegral o0) * ft0 (Z :. tf0) +
                 (deltaScale * fromIntegral s) * fs (Z :. sf) +
                 (deltaScale * fromIntegral s0) * fs0 (Z :. sf0))))
      mat2 =
        ((numOrientations * numScales * numOrientations * numScales) ><
         (xLen * yLen)) .
        R.toList $
        arr
   in fromListUnboxed
        (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :.
         xLen :.
         yLen) .
      NL.toList . flatten $
      mat1 NL.<> mat2
