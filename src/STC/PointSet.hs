module STC.PointSet where

import           Control.Monad as M
import           Data.List     as L
import           System.Random
import           Types

generateRandomPointSet :: Int -> Int -> Int -> IO [R2S1RPPoint]
generateRandomPointSet n rows cols =
  M.replicateM n $ do
    x <- randomRIO (0, cols - 1) :: IO Int
    y <- randomRIO (0, rows - 1) :: IO Int
    return . R2S1RPPoint $ (x, y, 0, 0)
    
generateRandomPointSet' :: Int -> (Double,Double) -> IO [(Double,Double)]
generateRandomPointSet' n range =
  M.replicateM n $ do
    x <- randomRIO range :: IO Double
    y <- randomRIO range :: IO Double
    return (x, y)
