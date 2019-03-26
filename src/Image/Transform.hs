{-# LANGUAGE FlexibleContexts #-}
module Image.Transform where

import           Data.Array.Repa              as R
import           Data.Array.Repa.Stencil      as R
import           Data.Array.Repa.Stencil.Dim2 as R
import           Data.DList                   as DL
import           Data.List                    as L
import           Data.Vector.Unboxed          as VU
import           Numeric.LinearAlgebra        as NL

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
pad :: (Source s e, Shape sh) => [Int] -> e -> R.Array s sh e -> R.Array D sh e
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

{-# INLINE normalizeValueRange #-}
normalizeValueRange ::
     (Shape sh, Num e, Fractional e, Ord e, Unbox e)
  => (e, e)
  -> R.Array U sh e
  -> R.Array U sh e
normalizeValueRange (minVal, maxVal) arr
  | maxV == minV = arr
  | otherwise =
    computeUnboxedS .
    R.map (\x -> (x - minV) / (maxV - minV) * (maxVal - minVal) + minVal) $
    arr
  where
    minV = VU.minimum . toUnboxed $ arr
    maxV = VU.maximum . toUnboxed $ arr

{-# INLINE computeDerivativeS #-}
computeDerivativeS ::
     (R.Source r e, Num e, Fractional e) => R.Array r DIM2 e -> R.Array D DIM3 e
computeDerivativeS arr =
  L.foldl' R.append (R.extend (Z :. R.All :. R.All :. (1 :: Int)) arr) ds
  where
    xStencil =
      makeStencil2 3 3 $ \ix ->
        case ix of
          Z :. 0 :. (-1) -> Just (-1)
          Z :. 0 :. 1 -> Just 1
          _ -> Nothing
    yStencil =
      makeStencil2 3 3 $ \ix ->
        case ix of
          Z :. -1 :. 0 -> Just (-1)
          Z :. 1 :. 0 -> Just 1
          _ -> Nothing
    xyStencil =
      makeStencil2 3 3 $ \ix ->
        case ix of
          Z :. -1 :. -1 -> Just 1
          Z :. -1 :. 1 -> Just (-1)
          Z :. 1 :. -1 -> Just (-1)
          Z :. 1 :. 1 -> Just 1
          _ -> Nothing
    ds =
      L.map
        (\s ->
           extend (Z :. R.All :. R.All :. (1 :: Int)) . R.map (/ 2) $
           mapStencil2 BoundClamp s arr)
        [xStencil, yStencil, xyStencil]
        
{-# INLINE bicubicInterpolation #-}
bicubicInterpolation ::
     (R.Source r Double)
  => ((Int, Int) -> (Double, Double))
  -> (Int, Int)
  -> R.Array r DIM2 Double
  -> R.Array U DIM2 Double
bicubicInterpolation idxFunc (newCols, newRows) arr =
  let (Z :. cols :. rows) = extent arr
      derivative = computeDerivativeS arr
      derivative16 =
        R.traverse
          derivative
          (const (Z :. newCols :. newRows :. 4 :. (4 :: Int)))
          (\f (Z :. a :. b :. derivativeIdx :. pairIdx) ->
             let (x, y) = idxFunc (a, b)
                 (x', y') =
                   case pairIdx of
                     0 -> (floor x, floor y)
                     1 -> (floor x, ceiling y)
                     2 -> (ceiling x, floor y)
                     3 -> (ceiling x, ceiling y)
              in if (x' < 0) ||
                    (x' > (fromIntegral cols - 1)) ||
                    (y' < 0) || (y' > (fromIntegral rows - 1))
                   then 0
                   else f (Z :. x' :. y' :. derivativeIdx))
      mat = ((newRows * newCols) >< 16) . R.toList $ derivative16
      multiplicationResultArray =
        fromListUnboxed (Z :. newCols :. newRows :. 4 :. 4) .
        NL.toList . flatten $
        mat NL.<> maxtrixA
   in R.sumS . R.sumS $
      R.traverse
        multiplicationResultArray
        id
        (\f idx@(Z :. a :. b :. i :. j) ->
           let (x, y) = idxFunc (a, b)
               x' = x - fromIntegral (floor x :: Int)
               y' = y - fromIntegral (floor y :: Int)
            in f idx * (x' ^ i) * (y' ^ j))
  where
    maxtrixA =
      NL.fromColumns . L.map NL.fromList $
      [ [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0]
      , [0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0]
      , [-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0]
      , [9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1]
      , [-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1]
      , [2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0]
      , [0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0]
      , [-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1]
      , [4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1]
      ] 

{-# INLINE resize2D #-}
resize2D ::
     (Source s Double)
  => (Int, Int)
  -> (Double, Double)
  -> R.Array s DIM2 Double
  -> R.Array U DIM2 Double
resize2D newSize@(newNx, newNy) bound arr =
  normalizeValueRange bound $
  bicubicInterpolation
    (\(i, j) -> (ratioX * fromIntegral i, ratioY * fromIntegral j))
    newSize
    arr
  where
    (Z :. nx' :. ny') = extent arr
    ratioX = fromIntegral nx' / fromIntegral newNx
    ratioY = fromIntegral ny' / fromIntegral newNy

{-# INLINE resize25D #-}
resize25D ::
     (Source s Double)
  => (Int, Int)
  -> (Double, Double)
  -> R.Array s DIM3 Double
  -> R.Array U DIM3 Double
resize25D newSize@(newNx, newNy) bound arr =
  normalizeValueRange bound .
  fromUnboxed (Z :. nf' :. newNx :. newNy) .
  VU.concat .
  L.map
    (\i ->
       toUnboxed .
       bicubicInterpolation
         (\(i, j) -> (ratioX * fromIntegral i, ratioY * fromIntegral j))
         newSize .
       R.slice arr $
       (Z :. i :. R.All :. R.All)) $
  [0 .. nf' - 1]
  where
    (Z :. nf' :. nx' :. ny') = extent arr
    ratioX = fromIntegral nx' / fromIntegral newNx
    ratioY = fromIntegral ny' / fromIntegral newNy
