{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE StrictData       #-}
module FourierPinwheel.Array where

import           Data.Array.Accelerate as A
import qualified Data.Array.Repa       as R
import           Data.Complex
import           Data.List             as L
import           Data.Vector.Generic   as VG
import           Data.Vector.Storable  as VS
import           DFT.Plan
import           Filter.Utils
import           Image.IO
import           Utils.Array
import           Utils.Parallel

data FPArray vector = FPArray
  { getFPArrayNumXFreq     :: Int
  , getFPArrayNumYFreq     :: Int
  , getFPArrayNumRFreq     :: Int
  , getFPArrayNumThetaFreq :: Int
  , getFPArrayNumRhoFreq   :: Int
  , getFPArrayNumPhiFreq   :: Int
  , getFPArray             :: [vector]
  }

{-# INLINE toMatrixAcc #-}
toMatrixAcc ::
     (VG.Vector vector e, Elt e) => FPArray (vector e) -> Acc (A.Array DIM2 e)
toMatrixAcc (FPArray numXFreq numYFreq numRFreq numThetaFreq _ _ arr) =
  A.use .
  A.fromList (Z :. (numRFreq * numThetaFreq) :. (numXFreq * numYFreq)) .
  VG.toList . VG.concat $
  arr

{-# INLINE parMapFPArray #-}
parMapFPArray ::
     (NFData vector) => (vector -> vector) -> FPArray vector -> FPArray vector
parMapFPArray f (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs) =
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
  parMap rdeepseq f $
  vecs

{-# INLINE parZipWithFPArray #-}
parZipWithFPArray ::
     (NFData vector)
  => (vector -> vector -> vector)
  -> FPArray vector
  -> FPArray vector
  -> FPArray vector
parZipWithFPArray f (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs1) arr =
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
  parZipWith rdeepseq f vecs1 . getFPArray $
  arr


plotFPArray ::
     DFTPlan
  -> FilePath
  -> FPArray (VS.Vector (Complex Double))
  -> IO (R.Array R.U R.DIM4 (Complex Double)) 
plotFPArray plan filePath arr = do
  let planIDBackward =
        DFTPlanID
          IDFT1DG
          [ getFPArrayNumThetaFreq arr
          , getFPArrayNumXFreq arr
          , getFPArrayNumYFreq arr
          ]
          [1, 2]
  vecsR2 <- dftExecuteBatchP plan planIDBackward . getFPArray $ arr
  arrR2 <-
    R.computeUnboxedP .
    makeFilter2DInverse .
    R.fromUnboxed
      (R.Z R.:. (getFPArrayNumRFreq arr) R.:. (getFPArrayNumThetaFreq arr) R.:.
       (getFPArrayNumXFreq arr) R.:.
       (getFPArrayNumYFreq arr)) .
    VS.convert . VS.concat $
    vecsR2
  arr1 <- R.sumP . R.sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $ arrR2
  plotImageRepa filePath .
    ImageRepa 8 .
    R.fromUnboxed
      (R.Z R.:. (1 :: Int) R.:. (getFPArrayNumXFreq arr) R.:.
       (getFPArrayNumYFreq arr)) .
    VG.map sqrt . R.toUnboxed $arr1
  return arrR2
