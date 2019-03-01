{-# LANGUAGE DeriveGeneric, DeriveAnyClass #-}
module FokkerPlanck.Histogram where

import           Control.DeepSeq
import           Data.Array.Repa     as R
import           Data.Binary
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           GHC.Generics        (Generic)

data Histogram a = Histogram
  { histogramSize :: [Int]
  , histogramNum  :: !Int
  , histogramVec  :: !(Vector a)
  } deriving (Generic,NFData)

instance (Binary a, Unbox a) => Binary (Histogram a) where
  put (Histogram size n vec) = do
    put size
    put n
    put . VU.toList $ vec
  get = do
    size <- get
    n <- get
    xs <- get
    return $! Histogram size n . VU.fromList $ xs

{-# INLINE addHistogram #-}
addHistogram :: (Num a, Unbox a) => Histogram a -> Histogram a -> Histogram a
addHistogram (Histogram size1 n1 vec1) (Histogram size2 n2 vec2) =
  if size1 == size2
    then Histogram size1 (n1 + n2) (VU.zipWith (+) vec1 vec2)
    else error $
         "addHistogram: sizes are not equal.\n" L.++ show size1 L.++ " /= " L.++
         show size2

{-# INLINE getNormalizedHistogramVec #-}
getNormalizedHistogramVec :: (Unbox a, Fractional a) => Histogram a -> VU.Vector a
getNormalizedHistogramVec (Histogram size n vec)
  | n == 0 = error "getNormalizedHistogramVec: n == 0."
  | (L.product size) /= (VU.length vec) =
    error $
    "getNormalizedHistogramVec: size mismatch\nSize = " L.++ show size L.++
    " = " L.++
    (show . L.product $ size) L.++
    "\nVector length = " L.++
    (show . VU.length $ vec)
  | otherwise = VU.map (/ fromIntegral n) vec

{-# INLINE getNormalizedHistogramArr #-}
getNormalizedHistogramArr ::
     (Unbox a, Fractional a, Shape d) => Histogram a -> R.Array U d a
getNormalizedHistogramArr (Histogram size n vec)
  | n == 0 = error "getNormalizedHistogramArr: n == 0."
  | (L.product size) /= (VU.length vec) =
    error $
    "getNormalizedHistogramArr: size mismatch\nSize = " L.++ show size L.++
    " = " L.++
    (show . L.product $ size) L.++
    "\nVector length = " L.++
    (show . VU.length $ vec)
  | otherwise = fromUnboxed (shapeOfList size) . VU.map (/ fromIntegral n) $ vec

{-# INLINE emptyHistogram #-}
emptyHistogram :: (Unbox a) => [Int] -> a -> Histogram a
emptyHistogram size = Histogram size 0 . VU.replicate (L.product size)

{-# INLINE mapHistogram #-}
mapHistogram :: (Unbox a, Unbox b) => (a -> b) -> Histogram a -> Histogram b
mapHistogram f (Histogram size n vec) = Histogram size n (VU.map f vec)
