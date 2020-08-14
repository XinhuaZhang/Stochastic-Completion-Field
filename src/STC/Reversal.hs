{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict #-}
module STC.Reversal where

import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           Data.Vector.Unboxed  as VU
import           STC.DFTArray
import           Types
import FourierPinwheel.Array

{-# INLINE computeReversalR2Z2T0S0 #-}
computeReversalR2Z2T0S0 ::
     (R.Source r (Complex Double))
  => [Double]
  -> R.Array r DIM6 (Complex Double)
  -> R.Array D DIM6 (Complex Double)
computeReversalR2Z2T0S0 theta0Freqs sourceArr =
  R.traverse2
    sourceArr
    (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
    const $ \f1 f2 idx@(Z :. _ :. _ :. l :. _ :. i :. j) ->
    f1 idx * exp (0 :+ f2 (Z :. l) * pi)

{-# INLINE timeReversal #-}
timeReversal :: DFTArray -> DFTArray
timeReversal arr@(DFTArray _ _ thetaFreqs rFreqs _) =
  repaToDFTArray thetaFreqs rFreqs .
  computeS .
  R.traverse2
    (dftArrayToRepa arr)
    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    const $ \fArr fFreq idx@(Z :. _ :. theta :. _ :. _) ->
    fArr idx * (cis $ (fFreq (Z :. theta)) * pi)

{-# INLINE timeReversal' #-}
timeReversal' :: DFTArray -> DFTArray
timeReversal' arr@(DFTArray rows cols thetaFreqs _ _) =
  repaToDFTArray thetaFreqs [] .
  computeS .
  R.traverse2
    (fromUnboxed (Z :. (1 :: Int) :. (L.length thetaFreqs) :. cols :. rows) .
     VS.convert . VS.concat . getDFTArrayVector $
     arr)
    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    const $ \fArr fFreq idx@(Z :. _ :. theta :. _ :. _) ->
    fArr idx * (cis $ (fFreq (Z :. theta)) * pi)
  
{-# INLINE timeReversal1 #-}
timeReversal1 :: Int -> Int -> DFTArray -> DFTArray
timeReversal1 numThetaFreq numRFreq arr@(DFTArray _ _ thetaFreqs rFreqs _) =
  repaToDFTArray thetaFreqs rFreqs .
  computeS .
  R.traverse2
    (dftArrayToRepa1 numThetaFreq numRFreq arr)
    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    const $ \fArr fFreq idx@(Z :. _ :. theta :. _ :. _) ->
    fArr idx * (cis $ (fFreq (Z :. theta)) * pi)

{-# INLINE timeReversalRepa #-}
timeReversalRepa ::
     (Num e, RealFloat e, R.Source s (Complex e), Unbox e)
  => [e]
  -> R.Array s DIM4 (Complex e)
  -> R.Array D DIM4 (Complex e)
timeReversalRepa thetaFreqs arr =
  R.traverse2
    arr
    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    const $ \fArr fFreq idx@(Z :. _ :. theta :. _ :. _) ->
    fArr idx * (cis $ (fFreq (Z :. theta)) * pi) 
