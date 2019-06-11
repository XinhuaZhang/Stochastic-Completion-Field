{-# LANGUAGE FlexibleContexts #-}
module STC.PowerMethodNormalization where

import           Data.Array.Repa      as R
import           Data.Complex
import           Data.Vector.Storable as VS
import           DFT.Plan
import           Filter.Utils
import           Types
import           Utils.Array

data PowerMethodNormalizationOption
  = PowerMethodPatchNorm (VS.Vector (Complex Double))
  | PowerMethodNone

{-# INLINE makePatchNormFilter #-}
makePatchNormFilter ::
     DFTPlan
  -> Int
  -> Int
  -> Bool
  -> Int
  -> IO (DFTPlan, PowerMethodNormalizationOption)
makePatchNormFilter plan cols rows True n = do
  let arr =
        makeFilter2D . fromFunction (Z :. cols :. rows) $ \(Z :. i :. j) ->
          if (sqrt . fromIntegral $ (i - div cols 2) ^ 2 + (j - div rows 2) ^ 2) <=
             fromIntegral n
            then 1
            else 0
  lock <- getFFTWLock
  (p, vec) <-
    dft1dGPlan lock plan [cols, rows] [0, 1] . VS.convert . toUnboxed . computeS $
    arr
  (newPlan, _) <- idft1dGPlan lock p [cols, rows] [0, 1] $ vec
  return (newPlan, PowerMethodPatchNorm vec)
makePatchNormFilter plan _ _ False _ = return (plan, PowerMethodNone)

{-# INLINE patchNorm #-}
patchNorm ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> VS.Vector (Complex Double)
  -> R.Array s DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
patchNorm plan filterF arr = do
  let (Z :. numThetaFreq :. numScaleFreq :. cols :. rows) = extent arr
  arrSum <- (R.sumP . rotate4D . rotate4D . R.map magnitude $ arr) >>= R.sumP
  arrSumF <-
    dftExecute plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    VS.convert . toUnboxed . computeS . R.map (:+ 0) $
    arrSum
  arrPatchSum <-
    fmap (fromUnboxed (Z :. cols :. rows) . VS.convert) .
    dftExecute plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) $
    VS.zipWith (*) arrSumF filterF
  computeP . R.traverse2 arr arrPatchSum const $ \f fNorm idx@(Z :. _ :. _ :. i :. j) ->
    if fNorm (Z :. i :. j) == 0
      then 0
      else f idx / fNorm (Z :. i :. j)

{-# INLINE powerMethodNormalization #-}
powerMethodNormalization ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> PowerMethodNormalizationOption
  -> R.Array s DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
powerMethodNormalization plan (PowerMethodPatchNorm filterF) arr =
  patchNorm plan filterF arr
powerMethodNormalization _ PowerMethodNone arr = do
  s <- R.sumAllP . rotate4D . rotate4D . R.map magnitude $ arr
  return . computeS . R.map (/ (s :+ 0)) $ arr
