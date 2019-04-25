{-# LANGUAGE FlexibleContexts #-}
module STC.Convolution where

import           Control.Monad.Parallel as MP (bindM2, mapM)
import           Data.Array.Repa        as R
import           Data.Complex
import           Data.List              as L
import           Data.Vector.Storable   as VS
import           DFT.Plan
import           Filter.Utils
import           Types

{-# INLINE makeFilterR2Z1 #-}
makeFilterR2Z1 ::
     (R.Source s (Complex Double))
  => R.Array s DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
makeFilterR2Z1 arr =
  let (Z :. numFreqs :. _ :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k) -> (Z :. (makeFilterHelper numFreqs i) :. j :. k))
        arr

{-# INLINE convolveR2T0 #-}
convolveR2T0 :: DFTPlan -> R2Z1T0Array -> R2T0Array -> IO R2Z1T0Array
convolveR2T0 plan filterF initialDistribution = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent filterF
  initialDistributionF <-
    fmap (fromUnboxed (Z :. numTheta0Freqs :. xLen :. yLen) . VS.convert) .
    dftExecute plan (DFTPlanID DFT1DG [numTheta0Freqs, xLen, yLen] [1, 2]) .
    VS.convert . toUnboxed $
    initialDistribution
  arrF <-
    computeP $
    R.traverse2
      filterF
      initialDistributionF
      const
      (\f1 f2 idx@(Z :. tf :. t0f :. x :. y) -> f1 idx * f2 (Z :. t0f :. x :. y))
  fmap (fromUnboxed (extent filterF) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreqs, numTheta0Freqs, xLen, yLen] [2, 3]) .
    VS.convert . toUnboxed $
    arrF


{-# INLINE makeFilterR2Z2 #-}
makeFilterR2Z2 ::
     (R.Source s (Complex Double))
  => R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
makeFilterR2Z2 arr =
  let (Z :. numThetaFreqs :. numScaleFreqs :. _ :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k :. l) ->
           (Z :. (makeFilterHelper numThetaFreqs i) :.
            (makeFilterHelper numScaleFreqs j) :.
            k :.
            l))
        arr

{-# INLINE convolveR2T0S0 #-}
convolveR2T0S0 :: DFTPlan -> R2Z2T0S0Array -> R2T0S0Array -> IO R2Z2T0S0Array
convolveR2T0S0 plan filterF initialDistribution = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent filterF
  initialDistributionF <-
    fmap
      (fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) .
       VS.convert) .
    dftExecute
      plan
      (DFTPlanID DFT1DG [numTheta0Freqs, numScale0Freqs, xLen, yLen] [2, 3]) .
    VS.convert . toUnboxed $
    initialDistribution
  arrF <-
    computeP $
    R.traverse2
      filterF
      initialDistributionF
      const
      (\f1 f2 idx@(Z :. tf :. sf :. t0f :. s0f :. x :. y) ->
         f1 idx * f2 (Z :. t0f :. s0f :. x :. y))
  fmap (fromUnboxed (extent filterF) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID
         IDFT1DG
         [ numThetaFreqs
         , numScaleFreqs
         , numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         ]
         [4, 5]) .
    VS.convert . toUnboxed $
    arrF


{-# INLINE convolveR2T0S0P #-}
convolveR2T0S0P ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> R2Z2T0S0Array
  -> R.Array s DIM4 (Complex Double)
  -> IO R2Z2T0S0Array
convolveR2T0S0P plan filterF initialDistribution = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent filterF
  initialDistributionF <-
    fmap
      (fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) .
       VS.convert . VS.concat) .
    MP.mapM
      (\t0f ->
         dftExecute plan (DFTPlanID DFT1DG [numScale0Freqs, xLen, yLen] [1, 2]) .
         VS.convert . toUnboxed . computeUnboxedS . R.slice initialDistribution $
         (Z :. t0f :. All :. All :. All)) $
    [0 .. numTheta0Freqs - 1]
  let arrF =
        R.traverse2
          filterF
          initialDistributionF
          const
          (\f1 f2 idx@(Z :. tf :. sf :. t0f :. s0f :. x :. y) ->
             f1 idx * f2 (Z :. t0f :. s0f :. x :. y))
  vecs <-
    MP.mapM
      (\(tf, sf) ->
         dftExecute
           plan
           (DFTPlanID
              IDFT1DG
              [numTheta0Freqs, numScale0Freqs, xLen, yLen]
              [2, 3]) .
         VS.convert . toUnboxed . computeS . R.slice arrF $
         (Z :. tf :. sf :. All :. All :. All :. All))
      [ (tf, sf)
      | tf <- [0 .. numThetaFreqs - 1]
      , sf <- [0 .. numScaleFreqs - 1]
      ]
  return . fromUnboxed (extent filterF) . VS.convert . VS.concat $ vecs

{-# INLINE timeReverseR2Z1 #-}
timeReverseR2Z1 ::
     (R.Source s (Complex Double))
  => [Double]
  -> R.Array s DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
timeReverseR2Z1 freqs arr =
  R.traverse2
    arr
    (fromListUnboxed (Z :. (L.length freqs)) freqs)
    const
    (\f1 f2 idx@(Z :. k :. i :. j) -> f1 idx * (exp (0 :+ f2 (Z :. k) * pi)))

{-# INLINE timeReversalConvolveR2Z1 #-}
timeReversalConvolveR2Z1 ::
     DFTPlan
  -> [Double]
  -> R.Array U DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
  -> IO (R.Array U DIM3 (Complex Double))
timeReversalConvolveR2Z1 plan thetaFreqs arr1 arr2 = do
  let (Z :. numThetaFreqs :. xLen :. yLen) = extent arr1
  vec1F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VS.convert . toUnboxed $
    arr1
  vec2F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VS.convert . toUnboxed . computeS . timeReverseR2Z1 thetaFreqs $
    arr2
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute plan (DFTPlanID IDFT1DG [numThetaFreqs, xLen, yLen] [0]) $
    VS.zipWith (*) vec1F vec2F

{-# INLINE convolveR2Z1 #-}
convolveR2Z1 ::
     DFTPlan
  -> [Double]
  -> R.Array U DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
  -> IO (R.Array U DIM3 (Complex Double))
convolveR2Z1 plan thetaFreqs arr1 arr2 = do
  let (Z :. numThetaFreqs :. xLen :. yLen) = extent arr1
  vec1F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VS.convert . toUnboxed $
    arr1
  vec2F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VS.convert . toUnboxed . computeS . makeFilterR2Z1 $
    arr2
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute plan (DFTPlanID IDFT1DG [numThetaFreqs, xLen, yLen] [0]) $
    VS.zipWith (\x y -> x * y) vec1F vec2F


{-# INLINE timeReverseR2Z2 #-}
timeReverseR2Z2 ::
     (R.Source s (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
timeReverseR2Z2 thetaFreqs scaleFreqs arr =
  R.traverse2
    arr
    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    const
    (\f1 f2 idx@(Z :. k :. l :. _ :. _) ->
       f1 idx * (exp (0 :+ f2 (Z :. k) * pi)))

{-# INLINE timeReversalConvolveR2Z2 #-}
timeReversalConvolveR2Z2 ::
     DFTPlan
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
timeReversalConvolveR2Z2 plan thetaFreqs scaleFreqs arr1 arr2 = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen) = extent arr1
  vec1F <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
    VS.convert . toUnboxed $
    arr1
  vec2F <-
    (computeP . makeFilterR2Z2 . timeReverseR2Z2 thetaFreqs scaleFreqs $ arr2) >>=
    (dftExecute 
       plan
       (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
     VS.convert . toUnboxed)
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) $
    VS.zipWith (*) vec1F vec2F


{-# INLINE convolveR2Z2 #-}
convolveR2Z2 ::
     DFTPlan
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
convolveR2Z2 plan thetaFreqs scaleFreqs arr1 arr2 = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen) = extent arr1
  vec1F <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
    VS.convert . toUnboxed $
    arr1
  vec2F <-
    (computeP . makeFilterR2Z2 $ arr2) >>=
    (dftExecute
       plan
       (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
     VS.convert . toUnboxed)
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) $
    VS.zipWith (*) vec1F vec2F
