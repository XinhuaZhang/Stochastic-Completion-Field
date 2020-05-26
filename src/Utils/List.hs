{-# LANGUAGE Strict #-}
module Utils.List where

import           Data.List as L

{-# INLINE divideList #-}
divideList :: Int -> [a] -> [[a]]
divideList n xs =
  if n <= 0
    then error $ "divideList: " L.++ (show n) L.++ " <= 0."
    else case L.splitAt n xs of
           (a, []) -> [a]
           (a, b) -> a : divideList n b
