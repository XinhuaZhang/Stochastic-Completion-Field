{-# LANGUAGE FlexibleContexts #-}
module STC.Utils where

import           Data.Array.Repa     as R
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
