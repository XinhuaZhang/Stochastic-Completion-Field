{-# LANGUAGE Strict #-}
module FourierPinwheel.Multiplication where

import           Control.Concurrent.Async
import           Data.Array.Repa          as R
import           Data.Complex
import           Data.List                as L
import           Data.Vector.Storable     as VS
import           DFT.Plan
import           FourierPinwheel.Array
import           Utils.Parallel
import Filter.Utils
import Utils.Distribution

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
    fromUnboxed (Z :. numRFreq :. numThetaFreq :. numXFreq :. numYFreq) .
    VS.convert <$>
    dftExecute plan planIDBackward vecs3F
  return .
    FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
    parMap
      rdeepseq
      (\rFreq ->
         VS.convert . toUnboxed . computeS . R.slice arr $ (Z :. rFreq :. All :. All :. All)) $
    [0 .. numRFreq - 1]


{-# INLINE multiplyBiasDiscrete #-}
multiplyBiasDiscrete ::
     DFTPlan
  -> VS.Vector (Complex Double)
  -> FPArray (VS.Vector (Complex Double))
  -> IO (FPArray (VS.Vector (Complex Double)))
multiplyBiasDiscrete plan bias (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs2F) = do
  let planIDForward = DFTPlanID DFT1DG [numThetaFreq, numXFreq, numYFreq] [1, 2]
      planIDBackward =
        DFTPlanID IDFT1DG [numThetaFreq, numXFreq, numYFreq] [1, 2]
  vecs2 <-
    parMap rdeepseq (VS.zipWith (*) bias) <$>
    dftExecuteBatchP plan planIDBackward vecs2F
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq <$>
    dftExecuteBatchP plan planIDForward vecs2
    

{-# INLINE multiplyRFunction #-}
multiplyRFunction ::
     DFTPlan
  -> Double
  -> FPArray (VS.Vector (Complex Double))
  -> IO (FPArray (VS.Vector (Complex Double)))
multiplyRFunction plan periodEnv (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs) = do
  let freqConst = 2 * pi / log periodEnv
      filterVec =
        computeUnboxedS . makeFilter1D . fromFunction (Z :. numRFreq) $ \(Z :. i') ->
          let i = fromIntegral (i' - div numRFreq 2)
           in 1 / (1 :+ freqConst * fromIntegral i) -- gaussian1D i 10 :+ 0 --1 / (1 :+ (freqConst * fromIntegral i)^2)
      vec = VS.convert . VS.concat $ vecs
  filterF <-
    dftExecute plan (DFTPlanID DFT1DG [numRFreq] [0]) . VS.convert . toUnboxed $
    filterVec
  vecF <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [numRFreq, numThetaFreq, numXFreq, numYFreq] [0])
      vec
  let arrF =
        fromUnboxed (Z :. numRFreq :. numThetaFreq :. numXFreq :. numYFreq) .
        VS.convert $
        vecF
      convolvedArrFR =
        R.traverse2
          arrF
          (fromUnboxed (Z :. numRFreq) . VS.convert $ filterF)
          const $ \fArr fVec idx@(Z :. r :. _ :. _ :. _) ->
          fArr idx * fVec (Z :. r)
  convolvedArrF <- computeUnboxedP convolvedArrFR
  convolvedVec <-
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numRFreq, numThetaFreq, numXFreq, numYFreq] [0]) .
    VS.convert . toUnboxed $
    convolvedArrF
  let convolvedArr =
        fromUnboxed (Z :. numRFreq :. numThetaFreq :. numXFreq :. numYFreq) .
        VS.convert $
        convolvedVec
      convolvedVecs =
        parMap
          rdeepseq
          (\i ->
             VS.convert . toUnboxed . computeS . R.slice convolvedArr $
             (Z :. i :. All :. All :. All))
          [0 .. numRFreq - 1]
  return
    (FPArray
       numXFreq
       numYFreq
       numRFreq
       numThetaFreq
       numRhoFreq
       numPhiFreq
       convolvedVecs) 

-- {-# INLINE multiplyRFunction #-}
-- multiplyRFunction :: FPArray (VS.Vector (Complex Double)) -> FPArray (VS.Vector (Complex Double))
-- multiplyRFunction (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs') =
--   let vecs = L.reverse vecs'
--       zeroVec = VS.replicate (numThetaFreq * numXFreq * numYFreq) 0
--       convolvedSGNPositive = L.init . L.scanl' (VS.zipWith (+)) zeroVec $ vecs
--       convolvedSGNNegative = L.tail . L.scanr (VS.zipWith (+)) zeroVec $ vecs
--       convolvedSGN =
--         parZipWith
--           rdeepseq
--           (VS.zipWith (\x y -> (x - y) * (0 :+ (-pi))))
--           convolvedSGNPositive
--           convolvedSGNNegative
--    in FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
--       L.reverse $
--       convolvedSGN
