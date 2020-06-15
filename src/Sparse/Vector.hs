{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE StrictData       #-}
module Sparse.Vector where

import           Control.DeepSeq
import           Data.Vector.Generic as VG

data SparseVector vector a = SparseVector
  { getSparseVectorIdx :: vector Int
  , getSparseVectorData :: vector a
  }
  
instance NFData (SparseVector vector a) where
  rnf (SparseVector vec1 vec2) = vec1 `seq` vec2 `seq` ()

{-# INLINE dot #-}
dot ::
     (VG.Vector vector e, Num e, VG.Vector vector Int)
  => SparseVector vector e
  -> vector e
  -> e
dot (SparseVector idx vec1) vec2 =
  VG.sum . VG.zipWith (*) vec1 $ VG.unsafeBackpermute vec2 idx


{-# INLINE getVectorIndex2D #-}
getVectorIndex2D :: Int -> Int -> Int -> Int
getVectorIndex2D cols i j = i * cols + j
