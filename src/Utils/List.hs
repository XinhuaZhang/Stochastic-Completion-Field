{-# LANGUAGE Strict #-}
module Utils.List where

import           Data.List as L

{-# INLINE divideList #-}
divideList :: Int -> Int -> [a] -> [[a]]
divideList 1 _ xs = [xs]
divideList n m xs =
  case L.splitAt m xs of
    (a, []) -> [a]
    (a, b) -> a : divideList (n - 1) m b

{-# INLINE divideListN #-}
divideListN :: Int -> [a] -> [[a]]
divideListN n xs =
  let len = fromIntegral . L.length $ xs :: Double
      m = ceiling $ len / (fromIntegral n) :: Int
  in if ceiling (len / fromIntegral m) < n
       then divideList n (m - 1) xs
       else divideList n m xs

{-# INLINE getListFromNumber #-}
getListFromNumber :: Int -> [Int]
getListFromNumber n =
  let m = div n 2
  in if odd n
       then [-m .. m]
       else [-m .. (m - 1)]
       
{-# INLINE getListFromNumber' #-}
getListFromNumber' :: (Floating a ) => Int -> [a]
getListFromNumber' n =
  let m = div n 2
   in if odd n
        then L.map fromIntegral [-m .. m]
        else L.map fromIntegral [-m .. (m - 1)]
