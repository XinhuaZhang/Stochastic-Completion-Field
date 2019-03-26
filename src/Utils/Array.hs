{-# LANGUAGE TypeOperators #-}
module Utils.Array where

import           Data.Array.Repa as R
import           Types

{-# INLINE rotate3D #-}
rotate3D :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e  
rotate3D arr =
  let (Z :. k :. i :. j) = extent arr
   in R.backpermute
        (Z :. i :. j :. k)
        (\(Z :. a :. b :. c) -> (Z :. c :. a :. b))
        arr
        

{-# INLINE rotate3DR #-}
rotate3DR :: (R.Source s e) => R.Array s DIM3 e -> R.Array D DIM3 e  
rotate3DR arr =
  let (Z :. i :. j :. k) = extent arr
   in R.backpermute
        (Z :. k :. i :. j)
        (\(Z :. a :. b :. c) -> (Z :. b :. c :. a))
        arr

{-# INLINE rotate4D #-}
rotate4D :: (R.Source s e) => R.Array s DIM4 e -> R.Array D DIM4 e  
rotate4D arr =
  let (Z :. k :. l :. i :. j) = extent arr
   in R.backpermute
        (Z :. l :. i :. j :. k)
        (\(Z :. a :. b :. c :. d) -> (Z :. d :. a :. b :. c))
        arr

{-# INLINE rotate4D' #-}
rotate4D' :: (R.Source s e) => R.Array s DIM4 e -> R.Array D DIM4 e  
rotate4D' arr =
  let (Z :. k :. l :. i :. j) = extent arr
   in R.backpermute
        (Z :. j :. k :. l :. i)
        (\(Z :. a :. b :. c :. d) -> (Z :. b :. c :. d :. a))
        arr

{-# INLINE rotateR2Z1T0Array #-}
rotateR2Z1T0Array :: (R.Source s e) => R.Array s DIM4 e -> R.Array D DIM4 e
rotateR2Z1T0Array arr =
  let (Z :. k :. l :. i :. j) = extent arr
   in R.backpermute
        (Z :. k :. i :. j :. l)
        (\(Z :. a :. b :. c :. d) -> (Z :. a :. d :. b :. c))
        arr
        
        

{-# INLINE rotateR2Z2T0S0Array #-}
rotateR2Z2T0S0Array :: (R.Source s e) => R.Array s DIM6 e -> R.Array D DIM6 e
rotateR2Z2T0S0Array arr =
  let (Z :. k :. l :. m :. n :. i :. j) = extent arr
  in R.backpermute
       (Z :. k :. l :. i :. j :. m :. n)
       (\(Z :. a :. b :. c :. d :. e :. f) -> (Z :. a :. b :. e :. f :. c :. d))
       arr
