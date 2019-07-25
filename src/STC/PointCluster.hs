module STC.PointCluster
  ( connectionMatrix
  , connectionMatrixP
  , createIndex2D
  , pointCluster
  , cluster2Array
  )
where

import           Data.Array.Repa    as R
import qualified Data.Array.Unboxed as AU
import           Data.List          as L
import           Utils.Parallel

type Index2D = (Int, (Int, Int))
type BoolArray2D = AU.Array (Int, Int) Bool

{-# INLINE isConnected #-}
isConnected :: Int -> (Int, Int) -> (Int, Int) -> Bool
isConnected 1 (x1, y1) (x2, y2) = abs (x1 - x2) <= 1 && abs (y1 - y2) <= 1
isConnected n (x1, y1) (x2, y2) =
  let r = sqrt . fromIntegral $ (x1 - x2) ^ 2 + (y1 - y2) ^ 2
   in if r <= fromIntegral n
        then True
        else False

{-# INLINE createIndex2D #-}
createIndex2D :: [(Int, Int)] -> [Index2D]
createIndex2D = L.zip [0..]

{-# INLINE connectionMatrix #-}
connectionMatrix :: Int -> [Index2D] -> BoolArray2D
connectionMatrix n xs =
  AU.array
      ((0, 0), (L.length xs - 1, L.length xs - 1))
      [((fst x, fst y), isConnected n (snd x) (snd y)) | x <- xs, y <- xs]

{-# INLINE connectionMatrixP #-}
connectionMatrixP ::
     ParallelParams -> Int -> [Index2D] -> BoolArray2D
connectionMatrixP parallelParams n xs =
  AU.array ((0, 0), (L.length xs - 1, L.length xs - 1)) .
  parMapChunk
    parallelParams
    rdeepseq
    (\(x, y) -> ((fst x, fst y), isConnected n (snd x) (snd y))) $
  [(x, y) | x <- xs, y <- xs]

{-# INLINE lookupMatrix #-}
lookupMatrix :: BoolArray2D -> Index2D -> Index2D -> Bool
lookupMatrix arr (i, _) (j, _) = arr AU.! (i, j)

{-# INLINE cluster1 #-}
cluster1 ::
     BoolArray2D
  -> Bool
  -> [Index2D]
  -> [Index2D]
  -> (Bool, [Index2D], [Index2D])
cluster1 _ flag xs [] = (flag, xs, [])
cluster1 arr flag xs (y:ys) =
  let (newFlag, connectedPoints, zs) = cluster1 arr flag xs ys
   in if L.any (lookupMatrix arr y) xs
        then cluster1 arr True (y : xs) ys
        else (newFlag, connectedPoints, y : zs)

{-# INLINE cluster #-}
cluster :: BoolArray2D -> [Index2D] -> [Index2D] -> ([Index2D], [Index2D])
cluster arr xs ys =
  let (flag, connectionedPoints, zs) = cluster1 arr False xs ys
   in if flag
        then cluster arr connectionedPoints zs
        else (connectionedPoints, zs)

pointCluster :: BoolArray2D -> [Index2D] -> [[(Int, Int)]]
pointCluster _ [] = []
pointCluster arr (x:xs) =
  let (ys, zs) = cluster arr [x] xs
   in L.map snd ys : (pointCluster arr zs)

{-# INLINE cluster2Array #-}
cluster2Array :: Int -> Int -> [[(Int, Int)]] -> [R.Array U DIM3 Double]
cluster2Array rows cols =
  L.map
    (\xs ->
       fromListUnboxed (Z :. (1 :: Int) :. cols :. rows) . AU.elems $
       (AU.accumArray (+) 0 ((0, 0), (cols - 1, rows - 1)) .
        L.map (\x -> (x, 255)) $
        xs :: AU.Array (Int, Int) Double))
