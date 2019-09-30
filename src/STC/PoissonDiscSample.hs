{-# LANGUAGE MonoLocalBinds #-}
{-# LANGUAGE FlexibleContexts #-}
module STC.PoissonDiscSample where

import           Control.Arrow
import           Control.Monad
import           Control.Monad.ST
import           Data.Array       as Arr
import           Data.Array.ST    as Arr
import           Data.Ix
import           Data.List        as L
import           Data.Sequence    as Seq
import           System.Random

data GridElement
  = GridElementIdx (Int, Int)
  | GridElementNone
  deriving (Show)

{-# INLINE gridElement2Index #-}
gridElement2Index :: [GridElement] -> [(Int,Int)]
gridElement2Index =
  L.map (\(GridElementIdx x) -> x) .
  L.filter
    (\x ->
       case x of
         GridElementIdx _ -> True
         GridElementNone  -> False)

{-# INLINE computeGridIndex #-}
computeGridIndex :: Double -> (Int, Int) -> (Int, Int)
computeGridIndex cellSize =
  join (***) (\x -> floor $ (fromIntegral x :: Double) / cellSize)

{-# INLINE computeGridIndexList #-}
computeGridIndexList :: Double -> ((Int,Int),(Int,Int)) -> (Int,Int) -> [(Int,Int)]
computeGridIndexList cellSize range' (x, y) =
  let (i, j) = computeGridIndex cellSize (x,y)
   in L.filter
        (inRange range')
        [(i + a, j + b) | a <- [-2 .. 2], b <- [-2 .. 2]]
        
searchCandidate ::
     (MArray (STArray s) GridElement (ST s), RandomGen g)
  => g
  -> Double
  -> ((Int, Int), (Int, Int))
  -> (Seq (Int, Int), STArray s (Int, Int) GridElement)
  -> ST s (Seq (Int, Int), (STArray s (Int, Int) GridElement))
searchCandidate gen cellSize range' (queue, arr) = do
  let (n, g) = randomR (0, Seq.length queue - 1) gen
      idx = computeGridIndexList cellSize range' . Seq.index queue $ n
  return undefined


-- sample :: Int -> Int -> Int -> Double -> [(Int, Int)] -> IO [(Int, Int)]
-- sample width height k radius (x:xs) =
--   let cellSize = radius / (sqrt 2)
--       gridWidth = ceiling $ fromIntegral width / cellSize
--       gridHeight = ceiling $ fromIntegral height / cellSize
--    in gridElement2Index . elems $
--       (runSTArray
--          (do grid <-
--                newArray
--                  ((0, 0), (gridHeight - 1, gridWidth - 1))
--                  GridElementNone
--              let queue = singleton x
--                  y = computeGridIndex cellSize x
--              writeArray grid y . GridElementIdx $ x
--              return grid) :: (Array (Int, Int) GridElement))
