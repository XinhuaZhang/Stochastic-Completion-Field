{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
module Utils.Diagram where

import           Control.Monad       as M
import           Data.Array.Repa     as R
import           Data.List           as L
import           Graphics.EasyRender
import           System.IO

{-# INLINE normalizeDirections #-}
normalizeDirections :: Double -> [[Double]] -> [[Double]]
normalizeDirections len xs =
  let !maxValue = L.maximum . L.concat $ xs
      !minValue = L.minimum . L.concat $ xs
      !y = len / (maxValue - minValue)
  in L.map (L.map (\x -> (x - minValue) * y)) xs
  
{-# INLINE normalizeDirections' #-}
normalizeDirections' :: Double -> [[Double]] -> [[Double]]
normalizeDirections' len =
  L.map
    (\xs ->
       let !maxValue = L.maximum $ xs
           !minValue = L.minimum $ xs
           !y = len / (maxValue - minValue)
       in L.map (\x -> (x - minValue) * y) xs)

drawR2S1 ::
     Double
  -> Double
  -> Double
  -> [(Double, Double)]
  -> [[Double]]
  -> Document ()
drawR2S1 !rows !cols !len xs directions =
  let normalizedDirections = normalizeDirections' len directions
      !rowCenter = rows / 2
      !colCenter = cols / 2
      !deltaTheat = 2 * pi / (fromIntegral . L.length . L.head $ directions)
  in newpage cols rows $
     M.zipWithM_
       (\(x, y) dirs -> do
          M.zipWithM_
            (\i v -> do
               let !theta = deltaTheat * fromIntegral i
                   !x' = v * cos theta
                   !y' = v * sin theta
               moveto (x + colCenter) (y + rowCenter)
               lineto (x + colCenter + x') (y + rowCenter + y')
               stroke)
            [0 ..]
            dirs)
       xs
       normalizedDirections

plotR2S1 ::
     FilePath
  -> Double
  -> Double
  -> Double
  -> [(Double, Double)]
  -> [[Double]]
  -> IO ()
plotR2S1 filePath !rows !cols !len xs directions =
  withFile filePath WriteMode $ \h ->
    render_file h (Format_EPS 1) (drawR2S1 rows cols len xs directions)

plotR2S1Array ::
     (R.Source r Double)
  => FilePath
  -> Double
  -> Double
  -> Double
  -> [(Double, Double)]
  -> R.Array r DIM3 Double
  -> IO ()
plotR2S1Array filePath !rows !cols !len xs arr =
  let !rowCenter = div (round rows) 2 :: Int
      !colCenter = div (round cols) 2 :: Int
      directions =
        L.map
          (\(x, y) ->
             R.toList . R.slice arr $
             (Z :. All :. (colCenter + round x) :. (rowCenter + round y)))
          xs
  in plotR2S1 filePath rows cols len xs directions
