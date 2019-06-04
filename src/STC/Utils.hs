{-# LANGUAGE FlexibleContexts #-}
module STC.Utils where

import           Control.Monad       as M
import           Data.Array.Repa     as R
import           Data.Array.ST       as Arr
import           Data.Array.Unboxed  as Arr hiding (Array)
import           Data.Ix
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           System.Random

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

{-# INLINE increasePixelDistance #-}
increasePixelDistance ::
     R.Source s Double => Int -> Array s DIM2 Double -> Array U DIM2 Double
increasePixelDistance n arr =
  fromListUnboxed (extent arr) . elems $
  ((runSTUArray $ do
      arrST <- newListArray ((0, 0), (rows - 1, cols - 1)) . R.toList $ arr
      M.mapM_
        (\(i, j) -> do
           x <- readArray arrST (i, j)
           when
             (x > 0)
             (do let idxs =
                       L.filter
                         (\(a, b) ->
                            inRange (0, rows - 1) a &&
                            inRange (0, cols - 1) b && (a, b) /= (i, j)) .
                       L.map (\(a, b) -> (i + a, j + b)) .
                       L.filter
                         (\(a, b) ->
                            (sqrt . fromIntegral $ a ^ 2 + b ^ 2) <
                            fromIntegral n) $
                       [(a, b) | a <- range, b <- range]
                 M.mapM_ (\idx -> writeArray arrST idx 0) idxs))
        [(i, j) | j <- [0 .. rows - 1], i <- [0 .. cols - 1]]
      return arrST) :: UArray (Int, Int) Double)
  where
    (Z :. rows :. cols) = extent arr
    range = [-n .. n]
    

{-# INLINE increasePixelDistanceRandom #-}
increasePixelDistanceRandom ::
     R.Source s Double => Int -> Array s DIM2 Double -> IO (Array U DIM2 Double)
increasePixelDistanceRandom n arr = do
  ys <- M.replicateM (rows * cols) randomIO :: IO [Double]
  let index =
        snd . L.unzip . L.sortOn fst . L.zip ys $
        [(i, j) | i <- [0 .. rows - 1], j <- [0 .. cols - 1]]
  return . fromListUnboxed (extent arr) . elems $
    ((runSTUArray $ do
        arrST <- newListArray ((0, 0), (rows - 1, cols - 1)) . R.toList $ arr
        M.mapM_
          (\(i, j) -> do
             x <- readArray arrST (i, j)
             when
               (x /= 0)
               (do let idxs =
                         L.filter
                           (\(a, b) ->
                              inRange (0, rows - 1) a &&
                              inRange (0, cols - 1) b && (a, b) /= (i, j)) .
                         L.map (\(a, b) -> (i + a, j + b)) .
                         L.filter
                           (\(a, b) ->
                              (sqrt . fromIntegral $ a ^ 2 + b ^ 2) <=
                              fromIntegral n) $
                         [(a, b) | a <- range, b <- range]
                   M.mapM_ (\idx -> writeArray arrST idx 0) idxs))
          index
        return arrST) :: UArray (Int, Int) Double)
  where
    (Z :. rows :. cols) = extent arr
    range = [-n .. n]
