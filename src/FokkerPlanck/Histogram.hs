{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric  #-}
module FokkerPlanck.Histogram where

import           Control.DeepSeq
import           Data.Array.Repa     as R
import           Data.Binary
import           Data.Complex
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           GHC.Generics        (Generic)
import           Types

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
  
{-# INLINE getHistogramArr #-}
getHistogramArr ::
     (Unbox a, Fractional a, Shape d) => Histogram a -> R.Array U d a
getHistogramArr (Histogram size n vec)
  | n == 0 = error "getHistogramArr: n == 0."
  | (L.product size) /= (VU.length vec) =
    error $
    "getHistogramArr: size mismatch\nSize = " L.++ show size L.++ " = " L.++
    (show . L.product $ size) L.++
    "\nVector length = " L.++
    (show . VU.length $ vec)
  | otherwise = fromUnboxed (shapeOfList size) vec

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
  
{-# INLINE getNormalizedHistogramArrSink #-}
getNormalizedHistogramArrSink ::
     (Shape d) => Histogram (Complex Double) -> R.Array U d (Complex Double)
getNormalizedHistogramArrSink (Histogram size n vec)
  | n == 0 = error "getNormalizedHistogramArr: n == 0."
  | (L.product size) /= (VU.length vec) =
    error $
    "getNormalizedHistogramArr: size mismatch\nSize = " L.++ show size L.++
    " = " L.++
    (show . L.product $ size) L.++
    "\nVector length = " L.++
    (show . VU.length $ vec)
  | otherwise =
    let normalizedVec =
          VU.map (\x -> (1 - (magnitude x) / fromIntegral n) :+ 0) $ vec
     in fromUnboxed (shapeOfList size) . VU.map (/ (VU.sum normalizedVec)) $
        normalizedVec
   

{-# INLINE getNormalizedHistogramArrSink' #-}
getNormalizedHistogramArrSink' ::
     (Shape d) => Histogram (Complex Double) -> R.Array U d (Complex Double)
getNormalizedHistogramArrSink' (Histogram size n vec)
  | n == 0 = error "getNormalizedHistogramArr: n == 0."
  | (L.product size) /= (VU.length vec) =
    error $
    "getNormalizedHistogramArr: size mismatch\nSize = " L.++ show size L.++
    " = " L.++
    (show . L.product $ size) L.++
    "\nVector length = " L.++
    (show . VU.length $ vec)
  | otherwise =
    let normalizedVec =
          VU.map
            (\x ->
               let m = magnitude x
                in if m == 0
                     then fromIntegral n :+ 0
                     else (fromIntegral n / m) :+ 0) $
          vec
     in fromUnboxed (shapeOfList size) -- . VU.map (/ (VU.sum normalizedVec)) $
        normalizedVec

{-# INLINE getNormalizedHistogramArrR2Z2T0S0 #-}
getNormalizedHistogramArrR2Z2T0S0 ::
     Histogram (Complex Double) -> R.Array U DIM5 (Complex Double)
getNormalizedHistogramArrR2Z2T0S0 (Histogram size n vec)
  | n == 0 = error "getNormalizedHistogramArr: n == 0."
  | (L.product size) /= (VU.length vec) =
    error $
    "getNormalizedHistogramArr: size mismatch\nSize = " L.++ show size L.++
    " = " L.++
    (show . L.product $ size) L.++
    "\nVector length = " L.++
    (show . VU.length $ vec)
  | otherwise =
    let arr = fromUnboxed (shapeOfList size) . VU.map (\x ->  magnitude x / fromIntegral n) $ vec
        maxArr = R.foldS max 0 arr
     in computeS . R.traverse2 arr maxArr const $ \f fMax idx@(Z :. a :. b :. c :. d :. r) ->
          (f (Z :. a :. b :. c :. d :. 1) - f idx) :+ 0


{-# INLINE emptyHistogram #-}
emptyHistogram :: (Unbox a) => [Int] -> a -> Histogram a
emptyHistogram size = Histogram size 0 . VU.replicate (L.product size)

{-# INLINE mapHistogram #-}
mapHistogram :: (Unbox a, Unbox b) => (a -> b) -> Histogram a -> Histogram b
mapHistogram f (Histogram size n vec) = Histogram size n (VU.map f vec)
