{-# LANGUAGE FlexibleContexts #-}
module STC.CompletionField where

import           Array.UnboxedArray   as AU
import           Control.Monad        as M
import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           Data.Vector.Unboxed  as VU
import           DFT.Plan
import           System.Random
import           Types

makeR2Z1T0Plan :: DFTPlan -> R2Z1T0Array -> IO DFTPlan
makeR2Z1T0Plan oldPlan (RepaArray arr) = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp2 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
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
                  ((x, y), exp (0 :+ (-1) * t0f * (theta / 180 * pi)))) $
             xs) $
        theta0Freqs
   in fmap
        (RepaArray .
         fromUnboxed (Z :. numTheta0Freqs :. xLen :. yLen) . VS.convert) .
      dftExecute plan (DFTPlanID DFT1DG [numTheta0Freqs, xLen, yLen] [1, 2]) .
      VS.convert $
      vec
      
{-# INLINE makeFilterR2Z1T0 #-}
makeFilterR2Z1T0 :: R2Z1T0Array -> R2Z1T0Array
makeFilterR2Z1T0 (RepaArray arr) =
  let (Z :. _ :. _ :. rows :. cols) = extent arr
   in RepaArray . computeS $
      R.backpermute
        (extent arr)
        (\(Z :. k :. l :. i :. j) ->
           let halfRows = div rows 2
               halfCols = div cols 2
               x =
                 if i < halfRows
                   then i + halfRows
                   else i - halfRows
               y =
                 if j < halfCols
                   then j + halfCols
                   else j - halfCols
            in (Z :. k :. l :. x :. y))
        arr
        
{-# INLINE makeFilterR2Z1 #-}
makeFilterR2Z1 ::
     (R.Source s (Complex Double))
  => R.Array s DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
makeFilterR2Z1 arr =
  let (Z :. numFreqs :. _ :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k) ->
           let halfNumFreqs = div numFreqs 2
               x =
                 if i < halfNumFreqs
                   then i + halfNumFreqs
                   else i - halfNumFreqs
            in (Z :. x :. j :. k))
        arr

{-# INLINE dftR2Z1T0 #-}
dftR2Z1T0 :: DFTPlan -> R2Z1T0Array -> IO R2Z1T0Array
dftR2Z1T0 plan (RepaArray arr) = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  fmap (RepaArray . fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numTheta0Freqs, xLen, yLen] [2, 3]) .
    VS.convert . toUnboxed $
    arr

{-# INLINE convolveR2T0 #-}
convolveR2T0 :: DFTPlan -> R2Z1T0Array -> R2T0Array -> IO R2Z1T0Array
convolveR2T0 plan (RepaArray filterF) (RepaArray initialDistributionF) = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent filterF
      arrF =
        computeS $
        R.traverse2
          filterF
          initialDistributionF
          const
          (\f1 f2 idx@(Z :. tf :. t0f :. x :. y) ->
             f1 idx * f2 (Z :. t0f :. x :. y))
  fmap (RepaArray . fromUnboxed (extent filterF) . VS.convert) .
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
    VU.convert .
    toUnboxed . computeS . makeFilterR2Z1 . timeReverseR2Z1 thetaFreqs $
    arr2
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute plan (DFTPlanID IDFT1DG [numThetaFreqs, xLen, yLen] [0]) $
    VS.zipWith (*) vec1F vec2F
