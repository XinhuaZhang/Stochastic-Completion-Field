{-# LANGUAGE TypeOperators #-}
module Filter.Utils where

import           Data.Array.Repa   as R
import           Data.List         as L


{-# INLINE center #-}
center :: Int -> Int
center x = div x 2

{-# INLINE makeFilterHelper #-}
makeFilterHelper :: Int -> Int -> Int
makeFilterHelper len i
  | even len =
    if i < c
      then i + c
      else i - c
  | otherwise =
    if i <= c
      then i + c
      else i - (c + 1)
  where
    c = center len

{-# INLINE makeFilter1D #-}
makeFilter1D ::
     (Source r e, Shape sh)
  => R.Array r (sh :. Int) e
  -> R.Array D (sh :. Int) e
makeFilter1D arr =
  let cols = L.head . listOfShape . extent $ arr
   in R.backpermute
        (extent arr)
        (\(sh :. i) -> (sh :. (makeFilterHelper cols i)))
        arr

{-# INLINE makeFilter2D #-}
makeFilter2D ::
     (Source r e, Shape sh)
  => R.Array r (sh :. Int :. Int) e
  -> R.Array D (sh :. Int :. Int) e
makeFilter2D arr =
  let (cols:rows:_) = L.take 2 . listOfShape . extent $ arr
   in R.backpermute
        (extent arr)
        (\(sh :. i :. j) ->
           (sh :. (makeFilterHelper rows i) :. (makeFilterHelper cols j)))
        arr
        
{-# INLINE makeFilter3D #-}
makeFilter3D ::
     (Source r e, Shape sh)
  => R.Array r (sh :. Int :. Int :. Int) e
  -> R.Array D (sh :. Int :. Int :. Int) e
makeFilter3D arr =
  let (sf:tf:cols:_) = L.take 3 . listOfShape . extent $ arr
   in R.backpermute
        (extent arr)
        (\(sh :. j :. a :. b) ->
           (sh :. (makeFilterHelper cols j) :. (makeFilterHelper tf a) :.
            (makeFilterHelper sf b)))
        arr
        

{-# INLINE makeFilter4D #-}
makeFilter4D ::
     (Source r e, Shape sh)
  => R.Array r (sh :. Int :. Int :. Int :. Int) e
  -> R.Array D (sh :. Int :. Int :. Int :. Int) e
makeFilter4D arr =
  let (sf:tf:cols:rows:_) = L.take 4 . listOfShape . extent $ arr
   in R.backpermute
        (extent arr)
        (\(sh :. i :. j :. a :. b) ->
           (sh :. (makeFilterHelper rows i) :. (makeFilterHelper cols j) :.
            (makeFilterHelper tf a) :.
            (makeFilterHelper sf b)))
        arr
