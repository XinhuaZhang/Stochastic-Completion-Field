{-# LANGUAGE FlexibleContexts #-}
module STC.PowerMethodNormalization where

import           Data.Array.Repa      as R
import           Data.Array.Unboxed   as AU
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           DFT.Plan
import           Filter.Utils
import           Types
import           Utils.Array
import           Utils.Parallel

data PowerMethodNormalizationOption
  = PowerMethodPatchNorm (VS.Vector (Complex Double))
  | PowerMethodGlobal
  | PowerMethodConnection [[(Int, Int)]]
  | PowerMethodConnectionMax [[(Int, Int)]]

-- {-# INLINE makePatchNormFilter #-}
-- makePatchNormFilter ::
--      DFTPlan
--   -> Int
--   -> Int
--   -> Bool
--   -> Int
--   -> IO (DFTPlan, PowerMethodNormalizationOption)
-- makePatchNormFilter plan cols rows True n = do
--   let arr =
--         makeFilter2D . fromFunction (Z :. cols :. rows) $ \(Z :. i :. j) ->
--           if (sqrt . fromIntegral $ (i - div cols 2) ^ 2 + (j - div rows 2) ^ 2) <=
--              fromIntegral n
--             then 1
--             else 0
--   lock <- getFFTWLock
--   (p, vec) <-
--     dft1dGPlan lock plan [cols, rows] [0, 1] . VS.convert . toUnboxed . computeS $
--     arr
--   (newPlan, _) <- idft1dGPlan lock p [cols, rows] [0, 1] $ vec
--   return (newPlan, PowerMethodPatchNorm vec)
-- makePatchNormFilter plan _ _ False _ = return (plan, PowerMethodGlobal)

-- {-# INLINE patchNorm #-}
-- patchNorm ::
--      (R.Source s (Complex Double))
--   => DFTPlan
--   -> VS.Vector (Complex Double)
--   -> R.Array s DIM4 (Complex Double)
--   -> IO (R.Array U DIM4 (Complex Double))
-- patchNorm plan filterF arr = do
--   let (Z :. numThetaFreq :. numScaleFreq :. cols :. rows) = extent arr
--   arrSum <- (R.sumP . rotate4D . rotate4D . R.map magnitude $ arr) >>= R.sumP
--   arrSumF <-
--     dftExecute plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
--     VS.convert . toUnboxed . computeS . R.map (:+ 0) $
--     arrSum
--   arrPatchSum <-
--     fmap (fromUnboxed (Z :. cols :. rows) . VS.convert) .
--     dftExecute plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) $
--     VS.zipWith (*) arrSumF filterF
--   computeP . R.traverse2 arr arrPatchSum const $ \f fNorm idx@(Z :. _ :. _ :. i :. j) ->
--     if fNorm (Z :. i :. j) == 0
--       then 0
--       else f idx / fNorm (Z :. i :. j)

-- {-# INLINE powerMethodNormalization #-}
-- powerMethodNormalization ::
--      (R.Source s (Complex Double))
--   => DFTPlan
--   -> PowerMethodNormalizationOption
--   -> R.Array s DIM4 (Complex Double)
--   -> IO (R.Array U DIM4 (Complex Double))
-- powerMethodNormalization plan (PowerMethodPatchNorm filterF) arr =
--   patchNorm plan filterF arr
-- powerMethodNormalization _ PowerMethodGlobal arr = do
--   s <- R.sumAllP . rotate4D . rotate4D . R.map magnitude $ arr
--   return . computeS . R.map (/ (s :+ 0)) $ arr

{-# INLINE powerMethodNormalization #-}
powerMethodNormalization ::
     (R.Source s (Complex Double))
  => PowerMethodNormalizationOption
  -> R.Array s DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
powerMethodNormalization PowerMethodGlobal arr = do
  -- s <- R.sumAllP . rotate4D . rotate4D . R.map (\x -> (magnitude x) ^ 2) $ arr
  s <- R.foldAllP max 0 . rotate4D . rotate4D . R.map (\x -> magnitude x) $ arr
  return . computeS . R.map (/ (s :+ 0)) $ arr
powerMethodNormalization (PowerMethodConnection xss) arr = do
  let (Z :. _ :. _ :. cols :. rows) = extent arr
      arrNorm =
        fromListUnboxed (Z :. cols :. rows) . AU.elems $
        (AU.accumArray (+) 0 ((0, 0), (cols - 1, rows - 1)) .
         L.concat .
         parMap
           rdeepseq
           (\xs ->
              let s =
                    L.maximum .
                    L.map
                      (\(i, j) ->
                         R.foldAllS max 0 .
                         R.slice (R.map (\x -> magnitude x) arr) $
                         (Z :. All :. All :. i :. j)) $
                    xs
               in L.map
                    (\x ->
                       ( x
                       , if s == 0
                           then 0
                           else 1 / s))
                    xs) $
         xss :: UArray (Int, Int) Double)
  return . computeS . R.traverse2 arr arrNorm const $ \f fNorm idx@(Z :. _ :. _ :. i :. j) ->
    f idx * ((fNorm (Z :. i :. j)) :+ 0)    
powerMethodNormalization (PowerMethodConnectionMax xss) arr = do
  let (Z :. numThetaFreq :. numScaleFreq :. cols :. rows) = extent arr
      arrMax =
        fromListUnboxed (Z :. cols :. rows) . AU.elems $
        (AU.accumArray (+) 0 ((0, 0), (cols - 1, rows - 1)) .
         L.concat .
         parMap
           rdeepseq
           (\xs ->
              let s =
                    sqrt .
                    L.maximum .
                    L.map
                      (\(i, j) ->
                         R.sumAllS .
                         R.slice (R.map (\x -> (magnitude x) ^ 2) arr) $
                         (Z :. All :. All :. i :. j)) $
                    xs
               in L.map (\x -> (x, s)) xs) $
         xss :: UArray (Int, Int) Double)
      arrSum =
        fromListUnboxed (Z :. cols :. rows) . AU.elems $
        (AU.accumArray (+) 0 ((0, 0), (cols - 1, rows - 1)) .
         L.concat .
         parMap
           rdeepseq
           (\xs ->
              L.map
                (\(i, j) ->
                   ( (i, j)
                   , sqrt .
                     R.sumAllS . R.slice (R.map (\x -> (magnitude x) ^ 2) arr) $
                     (Z :. All :. All :. i :. j))) $
              xs) $
         xss :: UArray (Int, Int) Double)
  return . computeS . R.traverse3 arr arrMax arrSum (\s1 s2 s3 -> s1) $ \f fMax fSum idx@(Z :. k :. l :. i :. j) ->
    if fSum (Z :. i :. j) > 0.9 * fMax (Z :. i :. j) && k == div numThetaFreq 2 && l == div numScaleFreq 2
      then 1 -- (f idx) / (fSum (Z :. i :. j) :+ 0)
      else 0

