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
  -> Double
  -> R.Array s DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
r2z2Tor2s1rp numOrientations thetaFreqs numScales scaleFreqs maxScale arr =
  let deltaTheta = 2 * pi / (fromIntegral numOrientations)
      deltaScale = log maxScale / (fromIntegral numScales)
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
                 (2 * pi * fromIntegral j * deltaScale / (log maxScale)) *
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
      deltaScale = log maxR / (fromIntegral numOrientations)
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
             (2 * pi * fromIntegral j * deltaScale / (log maxR)) * f3 (Z :. l))))
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

-- x y t s
{-# INLINE stsTor2s1rp #-}
stsTor2s1rp ::
     Int
  -> Double
  -> [Double]
  -> Int
  -> Double
  -> [Double]
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
stsTor2s1rp numX xLen xFreqs numY yLen yFreqs numOrientations thetaFreqs numScales scaleFreqs maxScale arr = do
  let deltaX = xLen / fromIntegral numX
      deltaY = yLen / fromIntegral numY
      deltaTheta = 2 * pi / fromIntegral numOrientations
      deltaScale = (log maxScale) / fromIntegral numScales
      centerX = div numX 2
      centerY = div numY 2
      vec = NL.fromList . R.toList $ arr
  transformMat <-
    fmap
      (((numX * numY * numOrientations * numScales) ><
        (L.product . listOfShape . extent $ arr)) .
       R.toList) .
    computeUnboxedP .
    R.traverse4
      (fromListUnboxed (Z :. L.length xFreqs) xFreqs)
      (fromListUnboxed (Z :. L.length yFreqs) yFreqs)
      (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
      (fromListUnboxed (Z :. L.length scaleFreqs) scaleFreqs)
      (\(Z :. numXFreq) (Z :. numYFreq) (Z :. numThetaFreq) (Z :. numScaleFreq) ->
         (Z :. numX :. numY :. numOrientations :. numScales :. numXFreq :.
          numYFreq :.
          numThetaFreq :.
          numScaleFreq)) $ \fx fy fo fs (Z :. x :. y :. o :. s :. xf :. yf :. orif :. sf) ->
      exp
        (0 :+
         (fo (Z :. orif) * fromIntegral o * deltaTheta +
          2 * pi *
          (fs (Z :. sf) * fromIntegral s * deltaScale / (log maxScale) +
           fx (Z :. xf) * fromIntegral ((x :: Int) - centerX) * deltaX / xLen +
           fy (Z :. yf) * fromIntegral ((y :: Int) - centerY) * deltaY / yLen)))
  return .
    fromListUnboxed (Z :. numX :. numY :. numOrientations :. numScales) .
    NL.toList $
    transformMat NL.#> vec
    

{-# INLINE stsTor2 #-}
stsTor2 ::
     Int
  -> Double
  -> [Double]
  -> Int
  -> Double
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM2 (Complex Double))
stsTor2 numX xLen xFreqs numY yLen yFreqs arr = do
  let deltaX = xLen / fromIntegral numX
      deltaY = yLen / fromIntegral numY
      centerX = div numX 2
      centerY = div numY 2
  vec <-
    fmap (NL.fromList . R.toList) .
    sumP . sumS . R.map (\x -> (magnitude x) ^ 2 :+ 0) $
    arr
  transformMat <-
    fmap
      (((numX * numY) >< (L.product . L.drop 2 . listOfShape . extent $ arr)) .
       R.toList) .
    computeUnboxedP .
    R.traverse2
      (fromListUnboxed (Z :. L.length xFreqs) xFreqs)
      (fromListUnboxed (Z :. L.length yFreqs) yFreqs)
      (\(Z :. numXFreq) (Z :. numYFreq) ->
         (Z :. numX :. numY :. numXFreq :. numYFreq)) $ \fx fy (Z :. x :. y :. xf :. yf) ->
      exp
        (0 :+
         (2 * pi *
          (fx (Z :. xf) * fromIntegral (x - centerX) * deltaX / xLen +
           fy (Z :. yf) * fromIntegral (y - centerY) * deltaY / yLen)))
  return . fromListUnboxed (Z :. numX :. numY) . NL.toList $
    transformMat NL.#> vec


{-# INLINE stsTor2' #-}
stsTor2' ::
     Int
  -> Double
  -> [Double]
  -> Int
  -> Double
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM2 Double)
stsTor2' numX xLen xFreqs numY yLen yFreqs arr = do
  let deltaX = xLen / fromIntegral numX
      deltaY = yLen / fromIntegral numY
      centerX = div numX 2
      centerY = div numY 2
      (Z :. numXF :. numYF :. numTF :. numSF) = extent arr
      mat = ((numXF * numYF) >< (numTF * numSF)) (R.toList arr)
  transformMat <-
    fmap
      (((numX * numY) >< (L.product . L.drop 2 . listOfShape . extent $ arr)) .
       R.toList) .
    computeUnboxedP .
    R.traverse2
      (fromListUnboxed (Z :. L.length xFreqs) xFreqs)
      (fromListUnboxed (Z :. L.length yFreqs) yFreqs)
      (\(Z :. numXFreq) (Z :. numYFreq) ->
         (Z :. numX :. numY :. numXFreq :. numYFreq)) $ \fx fy (Z :. x :. y :. xf :. yf) ->
      exp
        (0 :+
         (2 * pi *
          (fx (Z :. xf) * fromIntegral (x - centerX) * deltaX / xLen +
           fy (Z :. yf) * fromIntegral (y - centerY) * deltaY / yLen)))
  return .
    sumS .
    sumS .
    R.map (\x -> (magnitude x) ^ 2) .
    fromListUnboxed (Z :. numX :. numY :. numTF :. numSF) . NL.toList . flatten $
    transformMat NL.<> mat
