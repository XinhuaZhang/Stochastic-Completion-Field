{-# LANGUAGE Strict #-}
module FourierPinwheel.Multiplication where

import           Control.Concurrent.Async
import           Data.Array.Repa          as R
import           Data.Complex
import           Data.Vector.Storable     as VS
import           DFT.Plan
import           FourierPinwheel.Array
import           Utils.Parallel

{-# INLINE multiplyBias #-}
multiplyBias ::
     DFTPlan
  -> VS.Vector (Complex Double)
  -> FPArray (VS.Vector (Complex Double))
  -> IO (FPArray (VS.Vector (Complex Double)))
multiplyBias plan biasF (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs2) = do
  let planIDForward =
        DFTPlanID DFT1DG [numThetaFreq, numXFreq, numYFreq] [0, 1, 2]
      planIDBackward =
        DFTPlanID IDFT1DG [numThetaFreq, numXFreq, numYFreq] [0, 1, 2]
  vecs2F <- dftExecuteBatchP plan planIDForward vecs2
  let vecs3F = parMap rdeepseq (VS.zipWith (*) biasF) vecs2F
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq <$>
    dftExecuteBatchP plan planIDBackward vecs3F


{-# INLINE multiplyBias4D #-}
multiplyBias4D ::
     DFTPlan
  -> VS.Vector (Complex Double)
  -> FPArray (VS.Vector (Complex Double))
  -> IO (FPArray (VS.Vector (Complex Double)))
multiplyBias4D plan biasF (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs2) = do
  let planIDForward =
        DFTPlanID
          DFT1DG
          [numRFreq, numThetaFreq, numXFreq, numYFreq]
          [0, 1, 2, 3]
      planIDBackward =
        DFTPlanID
          IDFT1DG
          [numRFreq, numThetaFreq, numXFreq, numYFreq]
          [0, 1, 2, 3]
  vecs2F <- dftExecute plan planIDForward . VS.concat $ vecs2
  let vecs3F = VS.zipWith (*) biasF vecs2F
  arr <-
    fromUnboxed (Z :. numRFreq :. (numThetaFreq * numXFreq * numYFreq)) .
    VS.convert <$>
    dftExecute plan planIDBackward vecs3F
  return .
    FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
    parMap
      rdeepseq
      (\rFreq ->
         VS.convert . toUnboxed . computeS . R.slice arr $ (Z :. rFreq :. All)) $
    [0 .. numRFreq - 1]
