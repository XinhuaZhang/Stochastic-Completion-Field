{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}

module STC.CompletionField where

import           Array.UnboxedArray        as AU
import           Control.Monad             as M
import           Control.Monad.Parallel    as MP (bindM2,mapM)
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Storable      as VS
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.DomainChange (r2s1tor2z1, r2z1Tor2s1, r2z2Tor2s1rp, r2z2Tor2s1rpP)
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           Image.Transform           (normalizeValueRange)
import           System.FilePath           ((</>))
import           System.Random
import           Text.Printf
import           Types
import           Utils.Array
import           Utils.Time

makeR2Z1T0Plan ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> R.Array r DIM4 (Complex Double)
  -> IO DFTPlan
makeR2Z1T0Plan oldPlan arr = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp2 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp3 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [numThetaFreqs, numTheta0Freqs, xLen, yLen]
       [2, 3]
       vecTemp1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [numThetaFreqs, numTheta0Freqs, xLen, yLen]
         [2, 3]
         vec >>= \(plan, _) ->
         dft1dGPlan lock plan [numThetaFreqs, xLen, yLen] [0] vecTemp2 >>= \(plan, vec) ->
           idft1dGPlan lock plan [numThetaFreqs, xLen, yLen] [0] vec >>= \(plan, _) ->
             dft1dGPlan lock plan [numTheta0Freqs, xLen, yLen] [1, 2] vecTemp3 >>= \(plan, vec) ->
               idft1dGPlan lock plan [numTheta0Freqs, xLen, yLen] [1, 2] vec)


{-# INLINE dftR2Z1T0 #-}
dftR2Z1T0 :: DFTPlan -> R2Z1T0Array -> IO R2Z1T0Array
dftR2Z1T0 plan arr = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numTheta0Freqs, xLen, yLen] [2, 3]) .
    VS.convert . toUnboxed $
    arr

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
    VU.convert . toUnboxed $
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
    VU.convert . toUnboxed $
    arr1
  vec2F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VU.convert . toUnboxed . computeS . timeReverseR2Z1 thetaFreqs $
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
    VU.convert . toUnboxed $
    arr1
  vec2F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VU.convert . toUnboxed . computeS . makeFilterR2Z1 $
    arr2
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute plan (DFTPlanID IDFT1DG [numThetaFreqs, xLen, yLen] [0]) $
    VS.zipWith (\x y -> x * y) vec1F vec2F
    

{-# INLINE computeSinkFromSourceR2Z1T0 #-}
computeSinkFromSourceR2Z1T0 ::
     (R.Source r (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array r DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
computeSinkFromSourceR2Z1T0 thetaFreqs theta0Freqs sourceArr =
  R.traverse3
    sourceArr
    (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
    (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
    (\a _ _ -> a) $ \f1 f2 f3 idx@(Z :. k :. l :. i :. j) ->
    f1 idx * exp (0 :+ (f2 (Z :. k) + f3 (Z :. l)) * pi)

--- R2Z2T0S0 

makeR2Z2T0S0Plan ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> R.Array r DIM6 (Complex Double)
  -> IO DFTPlan
makeR2Z2T0S0Plan oldPlan arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent arr
  lock <- getFFTWLock
  let vecTemp1 = VU.convert . toUnboxed . computeS . delay $ arr
  fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [ numThetaFreqs
       , numScaleFreqs
       , numTheta0Freqs
       , numScale0Freqs
       , xLen
       , yLen
       ]
       [4, 5]
       vecTemp1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [ numThetaFreqs
         , numScaleFreqs
         , numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         ]
         [4, 5]
         vec >>= \(plan, vec) ->
         dft1dGPlan
           lock
           plan
           [numThetaFreqs, numScaleFreqs, xLen, yLen]
           [0, 1]
           vec >>= \(plan, vec) ->
           idft1dGPlan
             lock
             plan
             [numThetaFreqs, numScaleFreqs, xLen, yLen]
             [0, 1]
             vec >>= \(plan, vec) ->
             dft1dGPlan
               lock
               plan
               [numTheta0Freqs, numScale0Freqs, xLen, yLen]
               [2, 3]
               vec >>= \(plan, vec) ->
               idft1dGPlan
                 lock
                 plan
                 [numTheta0Freqs, numScale0Freqs, xLen, yLen]
                 [2, 3]
                 vec >>= \(plan, vec) ->
                 dft1dGPlan lock plan [numScale0Freqs, xLen, yLen] [1, 2] vec)

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
   in return $
      fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) vec
      
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
            yLen)) $ \f idx@(Z :. tf :. sf :. _ :. _ :. i :. j) ->
        if tf == div numThetaFreqs 2 && sf == div numScaleFreqs 2
          then f (Z :. i :. j) / fromIntegral (numTheta0Freqs * numScale0Freqs)
          else 0


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

{-# INLINE dftR2Z2T0S0 #-}
dftR2Z2T0S0 :: DFTPlan -> R2Z2T0S0Array -> IO R2Z2T0S0Array
dftR2Z2T0S0 plan arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent arr
  fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [ numThetaFreqs
         , numScaleFreqs
         , numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         ]
         [4, 5]) .
    VS.convert . toUnboxed $
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
         VU.convert . toUnboxed . computeUnboxedS . R.slice initialDistribution $
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
         VU.convert . toUnboxed . computeS . R.slice arrF $
         (Z :. tf :. sf :. All :. All :. All :. All))
      [ (tf, sf)
      | tf <- [0 .. numThetaFreqs - 1]
      , sf <- [0 .. numScaleFreqs - 1]
      ]
  return . fromUnboxed (extent filterF) . VS.convert . VS.concat $ vecs


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
    VU.convert . toUnboxed $
    arr1
  vec2F <-
    (computeP . makeFilterR2Z2 . timeReverseR2Z2 thetaFreqs scaleFreqs $ arr2) >>=
    (dftExecute 
       plan
       (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
     VU.convert . toUnboxed)
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
    VU.convert . toUnboxed $
    arr1
  vec2F <-
    (computeP . makeFilterR2Z2 $ arr2) >>=
    (dftExecute
       plan
       (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
     VU.convert . toUnboxed)
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) $
    VS.zipWith (*) vec1F vec2F
    
{-# INLINE computeSinkFromSourceR2Z2T0S0 #-}
computeSinkFromSourceR2Z2T0S0 ::
     (R.Source r (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array r DIM6 (Complex Double)
  -> R.Array D DIM6 (Complex Double)
computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceArr =
  R.traverse3
    sourceArr
    (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
    (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
    (\a _ _ -> a) $ \f1 f2 f3 idx@(Z :. k :. _ :. l :. _ :. i :. j) ->
    f1 idx * exp (0 :+ (f2 (Z :. k) + f3 (Z :. l)) * pi) 


-- Eigen methods
{-# INLINE eigenVectorR2Z1 #-}
eigenVectorR2Z1 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> R2Z1T0Array
  -> Int
  -> Bool
  -> String
  -> R.Array s1 DIM3 (Complex Double)
  -> R.Array s2 DIM4 (Complex Double)
  -> IO R2Z1T0Array
eigenVectorR2Z1 plan folderPath numOrientation thetaFreqs filterF n writeFlag name bias inputR2Z1T0 = do
  let (Z :. numThetaFreq :. numTheta0Freq :. cols :. rows) = extent inputR2Z1T0
  inputR2Z1 <- R.sumP . rotateR2Z1T0Array $ inputR2Z1T0
  let sourceDist = R.zipWith (*) bias inputR2Z1
      s = VU.maximum . VU.map magnitude . toUnboxed . computeS $ sourceDist
  let initialDist = computeS $ R.map (\x -> x / (s :+ 0)) sourceDist
  sourceArr <- convolveR2T0 plan filterF initialDist
  when
    writeFlag
    (do sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceArr
        sourceField <-
          computeP .
          R.extend (Z :. (1 :: Int) :. All :. All) .
          R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
          sourceR2Z1
        plotImageRepaComplex
          (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
          ImageRepa 8 $
          sourceField)
  return sourceArr

-- eigensink is computed from eigensource
powerMethod1 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> R2Z1T0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array s1 DIM3 (Complex Double)
  -> R.Array s2 DIM4 (Complex Double)
  -> IO (R.Array D DIM3 (Complex Double))
powerMethod1 plan folderPath cols rows numOrientation thetaFreqs theta0Freqs filter numIteration writeFlag idStr threshold bias eigenVecSource = do
  filterF <- dftR2Z1T0 plan . computeS . makeFilter2D $ filter
  sourceR2Z1T0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z1
           plan
           folderPath
           numOrientation
           thetaFreqs
           filterF
           n
           writeFlag
           "Source"
           bias
           input)
      (computeS . delay $ eigenVecSource)
      [1 .. numIteration]
  let sinkR2Z1T0 =
        computeSinkFromSourceR2Z1T0 thetaFreqs theta0Freqs sourceR2Z1T0
  sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceR2Z1T0
  sinkR2Z1 <- R.sumP . rotateR2Z1T0Array $ sinkR2Z1T0
  sinkField <-
    computeP .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sinkR2Z1
  plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  completionFieldR2 <-
    completionFieldR2Z1'
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      sourceR2Z1
      sinkR2Z1
  let completionFieldR2' =
        R.zipWith (\x y -> x * magnitude y) completionFieldR2 . R.slice bias $
        (Z :. (0 :: Int) :. All :. All)
      m = threshold * (VU.maximum . toUnboxed . computeS $ completionFieldR2)
      newBias =
        R.traverse2 bias completionFieldR2 const $ \fb fc idx@(Z :. _ :. i :. j) ->
          if fc (Z :. i :. j) >= m
            then 0
            else fb idx
  -- plotImageRepa
  --   (folderPath </> printf "bias%s.png" idStr)
  --   (ImageRepa 8 .
  --    computeS .
  --    R.map magnitude .
  --    R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
  --    newBias)
  return newBias

{-# INLINE timeReversal #-}

timeReversal :: (R.Source s e, Unbox e) => R.Array s DIM3 e -> R.Array U DIM3 e
timeReversal arr =
  let (Z :. nf :. _ :. _) = extent arr
      n = div nf 2
   in computeS $
      R.backpermute
        (extent arr)
        (\(Z :. k :. i :. j) ->
           let x = nf - n
            in if k >= x
                 then (Z :. k - x :. i :. j)
                 else (Z :. k + n :. i :. j))
        arr

completionFieldR2Z1 ::
     DFTPlan -> FilePath -> Int -> [Double]  -> R2T0Array -> R2T0Array -> IO ()
completionFieldR2Z1 plan folderPath numOrientation thetaFreqs source sink = do
  completionFiled <- timeReversalConvolveR2Z1  plan thetaFreqs (source) (sink)
  completionFiledR2 <-
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $ completionFiled
  plotImageRepaComplex (folderPath </> "TimeReversedCompletion.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> "TimeReversedCompletion_normalized.png") .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    reduceContrast 10 . R.map magnitude $
    completionFiledR2

completionFieldR2S1 ::
     DFTPlan
  -> FilePath
  -> String
  -> Int
  -> [Double]
  -> R2T0Array
  -> R2T0Array
  -> IO (R.Array D DIM2 Double)
completionFieldR2S1 plan folderPath idStr numOrientation thetaFreqs source sink = do
  let completionFiled =
        R.zipWith
          (*)
          (r2z1Tor2s1 numOrientation thetaFreqs $ source)
          (r2z1Tor2s1 numOrientation thetaFreqs $ sink)
  completionFiledR2 <- R.sumP . rotate3D $ completionFiled
  plotImageRepaComplex (folderPath </> printf "CompletioR2S1%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionR2S1_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map (logBase 10000) . normalizeValueRange (1, 256) . R.map magnitude $
    completionFiledR2
  let vec =
        toUnboxed .
        computeS .
        R.map (logBase 10000) . normalizeValueRange (1, 256) . R.map magnitude $
        completionFiledR2
  print . L.take 100 . L.reverse . L.sort . VU.toList $ vec
  return . R.map magnitude $ completionFiledR2

completionFieldR2Z1' ::
     DFTPlan -> FilePath -> String -> Int -> [Double]  -> R2T0Array -> R2T0Array -> IO (R.Array D DIM2 Double)
completionFieldR2Z1' plan folderPath idStr numOrientation thetaFreqs source sink = do
  completionFiled <- convolveR2Z1 plan thetaFreqs source sink
  completionFiledR2 <-
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $ completionFiled
  plotImageRepaComplex (folderPath </> printf "Completion%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionMagnitude%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
  let avg = R.sumAllS completionFiledR2 / (fromIntegral $ rows * cols)
      (Z :. cols :. rows) = extent completionFiledR2
  plotImageRepa (folderPath </> printf "Completion_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map log . normalizeValueRange (1, 256) . R.map magnitude $
    completionFiledR2
  return . R.map magnitude $ completionFiledR2
  

completionFieldR2Z2 ::
     DFTPlan
  -> FilePath
  -> String
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> R2T0S0Array
  -> R2T0S0Array
  -> IO (R.Array D DIM2 Double)
completionFieldR2Z2 plan folderPath idStr numOrientation thetaFreqs numScale scaleFreqs source sink = do
  completionFiled <- convolveR2Z2 plan thetaFreqs scaleFreqs source sink
  completionFiledR2 <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     completionFiled) >>=
    R.sumP
  plotImageRepaComplex (folderPath </> printf "Completion%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionMagnitude%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
  let avg = R.sumAllS completionFiledR2 / (fromIntegral $ rows * cols)
      (Z :. cols :. rows) = extent completionFiledR2
  plotImageRepa (folderPath </> printf "Completion_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map log . normalizeValueRange (1, 256) . R.map magnitude $
    completionFiledR2
  return . R.map magnitude $ completionFiledR2



{-# INLINE timeReverseR2Z1T0 #-}
timeReverseR2Z1T0 ::
     (R.Source s (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
timeReverseR2Z1T0 thetaFreqs theta0Freqs arr =
  let newArr =
        R.traverse3
          arr
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          (fromListUnboxed (Z :. (L.length theta0Freqs)) theta0Freqs)
          (\a _ _ -> a)
          (\f1 f2 f3 idx@(Z :. k :. l :. _ :. _) ->
             f1 idx * (exp (0 :+ (-f3 (Z :. l)) * pi)))
   in newArr

{-# INLINE timeReverseR2Z1T0' #-}
timeReverseR2Z1T0' ::
     (R.Source s (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
timeReverseR2Z1T0' thetaFreqs theta0Freqs arr =
  let newArr =
        R.traverse3
          arr
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          (fromListUnboxed (Z :. (L.length theta0Freqs)) theta0Freqs)
          (\a _ _ -> a)
          (\f1 f2 f3 idx@(Z :. k :. l :. _ :. _) ->
             f1 idx *
             (exp (0 :+ (f3 (Z :. l)  * pi / 2))))
      mag = R.sumAllS . R.map magnitude $ newArr
   in R.map (/ (mag :+ 0)) newArr




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

{-# INLINE eigenVectorR2Z2 #-}
eigenVectorR2Z2 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2 plan folderPath numOrientation thetaFreqs numScale scaleFreqs filterF n writeFlag name bias inputR2Z2T0S0 = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
        extent inputR2Z2T0S0
  inputR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ inputR2Z2T0S0) >>= R.sumP
  let sourceDist = R.zipWith (*) bias inputR2Z2
  s <- R.foldAllP max 0 . R.map magnitude $ sourceDist
  let initialDist = R.map (/ (s :+ 0)) sourceDist
  sourceArr <- convolveR2T0S0P plan filterF initialDist
  printCurrentTime $ printf "iteration %d" n
  when
    (n == 1 || (writeFlag && odd n))
    (let
      in do sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceArr) >>= R.sumP
            sourceR2S1RP <-
              r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
              sourceR2Z2
            sourceField <-
              (R.sumP . rotate4D . rotate4D $ sourceR2S1RP) >>=
              fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
              R.sumP
            plotImageRepaComplex
              (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
              ImageRepa 8 $
              sourceField)
  return sourceArr

powerMethodR2Z2T0S0 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> Int
  -> [Double]
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0 plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs filter numIteration writeFlag idStr threshold bias eigenVecSource = do
  filterF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filter
  sourceR2Z2T0S0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           filterF
           n
           writeFlag
           ("Source" L.++ idStr)
           bias
           input)
      (computeS . delay $ eigenVecSource)
      [1 .. numIteration]
  let sinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceR2Z2T0S0
  sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceR2Z2T0S0) >>= R.sumP
  sinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sinkR2Z2T0S0) >>= R.sumP
  sinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     sinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s.png" idStr) . ImageRepa 8 $
    sinkField
  completionFieldR2 <-
    completionFieldR2Z2
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      sourceR2Z2
      sinkR2Z2
  let m = threshold * (VU.maximum . toUnboxed . computeS $ completionFieldR2)
      newBias =
        R.traverse2 bias completionFieldR2 const $ \fb fc idx@(Z :. _ :. _ :. i :. j) ->
          if fc (Z :. i :. j) >= m
            then 0
            else fb idx
  -- plotImageRepa
  --   (folderPath </> printf "bias%s.png" idStr)
  --   (ImageRepa 8 .
  --    computeS .
  --    R.map magnitude .
  --    R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
  --    newBias)
  return newBias    

-- Utils

{-# INLINE dropPixel #-}
dropPixel :: Double -> VU.Vector Double -> IO (VU.Vector Double)
dropPixel t =
  VU.mapM
    (\x -> do
       v <- randomIO
       return $
         if v < t
           then 0
           else x)

{-# INLINE filterImage #-}
filterImage ::
     (R.Source s Double, Shape sh) => R.Array s sh Double -> R.Array D sh Double
filterImage arr =
  let avg = (R.sumAllS arr) / (fromIntegral . R.size . extent $ arr)
   in R.map
        (\x ->
           if x > avg
             then x
             else 0)
        arr
        
makeImagePlan ::
     DFTPlan
  -> R.Array U DIM3 Double
  -> IO (DFTPlan, R.Array U DIM3 (Complex Double))
makeImagePlan plan arr = do
  let (Z :. channels :. cols :. rows) = extent arr
  lock <- getFFTWLock
  (plan1, imgF) <-
    dft1dGPlan lock plan [channels, cols, rows] [1, 2] .
    VU.convert . VU.map (:+ 0) . toUnboxed $
    arr
  (newPlan, _) <- idft1dGPlan lock plan1 [channels, cols, rows] [1, 2] imgF
  return (newPlan, fromUnboxed (extent arr) . VS.convert $ imgF)


{-# INLINE reduceContrast #-}
reduceContrast ::
     (R.Source s Double, Shape sh)
  => Int
  -> Array s sh Double
  -> Array D sh Double
reduceContrast n arr =
  let x = L.head . L.drop n . L.reverse . L.sort . R.toList $ arr
   in R.map
        (\y ->
           if y >= x
             then x
             else y)
        arr
