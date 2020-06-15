{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Strict #-}
module Utils.SimpsonRule
  ( weightsSimpsonRule
  , computeWeightArrFromListOfShape
  , computeWeightArr
  , weightedArray
  , integralArray
  , integralArrayP
  ) where

import           Data.Array.Repa     as R
import           Data.List           as L
import           Data.Vector.Unboxed (Unbox)
import           Text.Printf

weightsHelper :: (Num a) => Int -> [a]
weightsHelper !0 = [1]
weightsHelper !n = 2 : (4 : weightsHelper (n - 1))

{-# INLINE weightsSimpsonRule #-}
weightsSimpsonRule :: (Num a) => Int -> [a]
weightsSimpsonRule n
  | (odd n) && (n >= 3) = 1 : (4 : (weightsHelper $ div (n - 3) 2))
  | otherwise =
    error $
    printf
      "weightsSimpsonRule: The number of weights (n = %d) for Simpson's rule ought to be odd and >= 3."
      n

{-# INLINE weightFunc #-}
weightFunc :: (Num a) => Int -> Int -> a
weightFunc !n !i
  | i == 0 || i == n - 1 = 1
  | odd i = 4
  | otherwise = 2
  
-- {-# INLINE weightFuncEven #-}
-- weightFuncEven :: (Num a) => Int -> Int -> a
-- weightFuncEven n i 
--   | i == 0 || i == n - 1 = 1
--   | rem == 0 = 2
--   | otherwise = 3
--   where rem = mod i 3


-- {-# INLINE weightFunc #-}
-- weightFunc :: (Num a) => Int -> Int -> a
-- weightFunc n i
--   | odd n = weightFuncOdd n i
--   | otherwise = weightFuncEven n i

{-# INLINE computeWeightArrFromListOfShape #-}
computeWeightArrFromListOfShape :: (Num a, Shape sh) => [Int] -> R.Array D sh a
computeWeightArrFromListOfShape xs =
  if L.any even xs || L.any (< 3) xs
    then error $
         printf
           "computeWeightArr: The number of weights for Simpson's rule ought to be odd and >= 3. %s" .
         show $
         xs
    else fromFunction (shapeOfList xs) $
         L.product . L.zipWith weightFunc xs . listOfShape 

{-# INLINE computeWeightArr #-}
computeWeightArr ::
     (Num a, R.Source r a, Shape sh) => R.Array r sh a -> R.Array D sh a
computeWeightArr = computeWeightArrFromListOfShape . listOfShape . extent

{-# INLINE weightedArray #-}
weightedArray ::
     (Num a, R.Source r a, Shape sh, Unbox a) => R.Array r sh a -> R.Array D sh a
weightedArray arr = arr *^ (computeUnboxedS $ computeWeightArr arr)

{-# INLINE integralArray #-}
integralArray ::
     (Num a, Fractional a, R.Source r a, Shape sh, Unbox a) => a -> R.Array r sh a -> a
integralArray delta arr =
  let !dim = L.length . listOfShape . extent $ arr
  in (delta / 3) ^ dim * sumAllS (weightedArray arr)

{-# INLINE integralArrayP #-}
integralArrayP ::
     (Num a, Fractional a, Unbox a, R.Source r a, Shape sh)
  => a
  -> R.Array r sh a
  -> IO a
integralArrayP delta arr =
  let !dim = L.length . listOfShape . extent $ arr
  in fmap (* (delta / 3) ^ dim) . sumAllP . weightedArray $ arr
