{-# LANGUAGE FlexibleContexts #-}
module STC.InitialDistribution where

import           Array.UnboxedArray  as AU
import           Control.Monad       as M
import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           DFT.Plan
import           STC.PointSet
import           System.Random
import           Types
import Utils.Coordinates

computeInitialDistributionR2T0 ::
     DFTPlan -> Int -> Int -> [Double] -> [R2S1RPPoint] -> IO R2T0Array
computeInitialDistributionR2T0 plan xLen yLen theta0Freqs xs =
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
        VU.concat .
        L.map
          (\t0f ->
             toUnboxedVector .
             AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
             L.map
               (\(R2S1RPPoint (x, y, theta, _)) ->
                  ( (x, y)
                  , (1 / (fromIntegral . L.length $ xs) :+ 0) *
                    exp (0 :+ (-1) * t0f * (theta / 180 * pi)))) $
             xs) $
        theta0Freqs
   in return $ fromUnboxed (Z :. numTheta0Freqs :. xLen :. yLen) vec

computeInitialEigenVectorR2T0 ::
     Int -> Int -> [Double] -> [Double] -> [R2S1RPPoint] -> R2Z1T0Array
computeInitialEigenVectorR2T0 xLen yLen theta0Freqs thetaFreqs xs =
  let numTheta0Freqs = L.length theta0Freqs
      numThetaFreqs = L.length thetaFreqs
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
        L.map
          (\(R2S1RPPoint (x, y, theta, _)) ->
             ((x, y), 1 / (fromIntegral . L.length $ xs))) $
        xs
   in computeS .
      R.traverse
        (fromUnboxed (Z :. xLen :. yLen) vec)
        (const (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen)) $ \f idx@(Z :. k :. _ :. i :. j) ->
        if k == div numThetaFreqs 2
          then f (Z :. i :. j) / fromIntegral numTheta0Freqs
          else 0


computeInitialDistributionR2T0S0 ::
     DFTPlan
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> Double
  -> [R2S1RPPoint]
  -> IO R2T0S0Array
computeInitialDistributionR2T0S0 plan xLen yLen theta0Freqs scale0Freqs maxScale xs =
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
               (\(R2S1RPPoint (x, y, theta, scale')) ->
                  let scale =
                        if scale' == 0
                          then 0
                          else log scale'
                   in ( (x, y)
                      , exp
                          (0 :+
                           (-1) *
                           (t0f * (theta / 180 * pi) +
                            s0f * 2 * pi * scale / log maxScale)))) $
             xs) $
        [(t0f, s0f) | t0f <- theta0Freqs, s0f <- scale0Freqs]
   in return . fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) $ vec
      -- computeS .
      -- R.traverse
      --   (fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) vec)
      --   id $ \f idx@(Z :. a :. b :. c :. d) ->
      --   if a == div (L.length theta0Freqs) 2 &&
      --      b == div (L.length scale0Freqs) 2
      --     then f idx
      --     else 0


computeInitialEigenVectorR2T0S0 ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> [R2S1RPPoint]
  -> R2Z2T0S0Array
computeInitialEigenVectorR2T0S0 xLen yLen theta0Freqs scale0Freqs thetaFreqs scaleFreqs xs =
  let numTheta0Freqs = L.length theta0Freqs
      numScale0Freqs = L.length scale0Freqs
      numThetaFreqs = L.length thetaFreqs
      numScaleFreqs = L.length scaleFreqs
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
        L.map
          (\(R2S1RPPoint (x, y, _, _)) ->
             ((x, y), 1 / (fromIntegral . L.length $ xs))) $
        xs
   in computeS .
      R.traverse
        (fromUnboxed (Z :. xLen :. yLen) vec)
        (const
           (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :.
            numScale0Freqs :.
            xLen :.
            yLen)) $ \f idx@(Z :. tf :. sf :. t0f :. s0f :. i :. j) -> 
        if tf == div numThetaFreqs 2 && sf == div numScaleFreqs 2
           -- t0f == div numThetaFreqs 2 && s0f == div numScaleFreqs 2
          then f (Z :. i :. j) / fromIntegral (numTheta0Freqs * numScale0Freqs)
          else 0
          

{-# INLINE computeInitialEigenVectorR2T0S0FromRepa #-}
computeInitialEigenVectorR2T0S0FromRepa ::
     (R.Source s Double)
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> R.Array s DIM2 Double
  -> R.Array D DIM6 (Complex Double)
computeInitialEigenVectorR2T0S0FromRepa rows cols numTheta0Freq numScale0Freq numThetaFreq numScaleFreq arr =
  R.traverse
    arr
    (const
       (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :.
        cols :.
        rows)) $ \f (Z :. k :. l :. _ :. _ :. i :. j) ->
    if k == div numThetaFreq 2 && l == div numScaleFreq 2 && f (Z :. i :. j) > 0
      then f (Z :. i :. j) / (fromIntegral $ (numTheta0Freq * numScale0Freq)) :+
           0
      else 0


{-# INLINE initializeEigenVectorR2Z1 #-}
initializeEigenVectorR2Z1 ::
     Int -> Int -> [Double] -> [Double] -> IO R2Z1T0Array
initializeEigenVectorR2Z1 xLen yLen thetaFreqs theta0Freqs = do
  let size = xLen * yLen
      numTheta0Freq = L.length theta0Freqs
  mags <- M.replicateM size randomIO
  thetas <- M.replicateM size (randomRIO (0, 2 * pi))
  return .
    computeS .
    R.traverse3
      (fromListUnboxed (Z :. xLen :. yLen) mags)
      (fromListUnboxed (Z :. xLen :. yLen) thetas)
      (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
      (\_ _ (Z :. numThetaFreq) ->
         (Z :. numThetaFreq :. numTheta0Freq :. xLen :. yLen)) $ \fMag fTheta fFreq idx@(Z :. k :. _ :. i :. j) ->
    (fMag (Z :. i :. j) / (fromIntegral numTheta0Freq) :+ 0) *
    exp (0 :+ (-1) * fFreq (Z :. k) * fTheta (Z :. i :. j))


-- For local eigenvector 
computeInitialEigenVector4D ::
     Int -> Int -> [Double] -> [Double] -> [R2S1RPPoint] -> R2Z2Array
computeInitialEigenVector4D xLen yLen thetaFreqs scaleFreqs xs =
  let numThetaFreqs = L.length thetaFreqs
      numScaleFreqs = L.length scaleFreqs
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
        L.map
          (\(R2S1RPPoint (x, y, _, _)) ->
             ((x, y), 1 / (fromIntegral . L.length $ xs))) $
        xs
   in computeS .
      R.traverse
        (fromUnboxed (Z :. xLen :. yLen) vec)
        (const (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen)) $ \f idx@(Z :. tf :. sf :. i :. j) ->
          f (Z :. i :. j)
        -- if tf == div numThetaFreqs 2 && sf == div numScaleFreqs 2
        --   then f (Z :. i :. j)
        --   else 0


computeInitialDistributionR2T0S0' ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> Double
  -> [R2S1RPPoint]
  -> R2T0S0Array
computeInitialDistributionR2T0S0' xLen yLen theta0Freqs scale0Freqs maxScale xs =
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
               (\(R2S1RPPoint (x, y, _, scale')) ->
                  let scale =
                        if scale' == 0
                          then 0
                          else log scale'
                      theta = angleFunctionRad (fromIntegral x) (fromIntegral y)
                   in ( (x, y)
                      , exp
                          (0 :+
                           (-1) *
                           (t0f * (0) +
                            s0f * 2 * pi * scale / log maxScale)) -- +
                        -- (exp
                        --    (0 :+
                        --     (-1) *
                        --     (t0f * (theta + pi ) +
                        --      s0f * 2 * pi * scale / log maxScale)))
                       )) $
             xs) $
        [(t0f, s0f) | t0f <- theta0Freqs, s0f <- scale0Freqs]
   in fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) $ vec
