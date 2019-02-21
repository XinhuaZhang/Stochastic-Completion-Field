module Image.Transform where

import           Data.Array.Repa     as R
import           Data.List           as L
import           Data.Vector.Unboxed as VU

-- factor = 2^n, n = 0,1,..
-- the first factor in the list corresponds to the inner-most (right-most) dimension.
{-# INLINE downsample #-}
downsample ::
     (Source s e, Shape sh) => [Int] -> R.Array s sh e -> R.Array D sh e
downsample factorList arr
  | L.all (== 1) factorList = delay arr
  | L.any (< 1) newDList = error "Downsample factors are too large."
  | otherwise =
    R.backpermute
      (shapeOfList newDList)
      (shapeOfList . L.zipWith (*) factorList . listOfShape)
      arr
  where
    dList = listOfShape . extent $ arr
    newDList = L.zipWith div dList factorList

{-# INLINE downsampleUnsafe #-}
downsampleUnsafe ::
     (Source s e, Shape sh) => [Int] -> R.Array s sh e -> R.Array D sh e
downsampleUnsafe factorList arr =
  R.backpermute newSh (shapeOfList . L.zipWith (*) factorList . listOfShape) arr
  where
    dList = listOfShape $ extent arr
    newSh = shapeOfList $ L.zipWith div dList factorList

{-# INLINE crop #-}
crop ::
     (Source s e, Shape sh)
  => [Int]
  -> [Int]
  -> R.Array s sh e
  -> R.Array D sh e
crop start len arr
  | L.any (< 0) start ||
      L.or (L.zipWith3 (\x y z -> x > (z - y)) start len dList) =
    error $
    "Crop out of boundary!\n" L.++ show start L.++ "\n" L.++ show len L.++ "\n" L.++
    show dList
  | L.length start /= L.length len || L.length start /= L.length dList =
    error $
    "crop: dimension error. \n start: " L.++ show (L.length start) L.++ " len:" L.++
    show (L.length len) L.++
    " arr:" L.++
    show (L.length dList)
  | otherwise =
    R.backpermute
      (shapeOfList len)
      (shapeOfList . L.zipWith (+) start . listOfShape)
      arr
  where
    dList = listOfShape $ extent arr

{-# INLINE cropUnsafe #-}
cropUnsafe ::
     (Source s e, Shape sh)
  => [Int]
  -> [Int]
  -> R.Array s sh e
  -> R.Array D sh e
cropUnsafe start len =
  R.backpermute
    (shapeOfList len)
    (shapeOfList . L.zipWith (+) start . listOfShape)

{-# INLINE pad #-}
pad ::
     (Real e, Source s e, Shape sh)
  => [Int]
  -> e
  -> R.Array s sh e
  -> R.Array D sh e
pad newDims padVal arr
  | L.all (== 0) diff = delay arr
  | otherwise =
    backpermuteDft
      (fromFunction (shapeOfList dimList) (const padVal))
      (\sh' ->
         let idx = L.zipWith (-) (listOfShape sh') diff
          in if L.or (L.zipWith (\i j -> i < 0 || (i >= j)) idx oldDimList)
               then Nothing
               else Just $ shapeOfList idx)
      arr
  where
    oldDimList = listOfShape . extent $ arr
    dimList = L.zipWith max newDims oldDimList
    diff =
      L.zipWith
        (\a b ->
           if a - b <= 0
             then 0
             else if a - b == 1
                    then 1
                    else div (a - b) 2)
        newDims
        oldDimList
