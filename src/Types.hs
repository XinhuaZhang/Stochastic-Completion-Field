module Types where

import           Data.Array.Repa
import           Data.Complex
import           Data.Vector.Storable as VS
import           Data.Vector.Unboxed  as VU

data DFTArray dim e
  = RepaArray (Array U dim e)
  | VectorArray [Int]
                (VS.Vector e)
                
type R2S1Array   = Array U DIM3 (Complex Double) -- (Z :. numOrientations :. xLen :. yLen)
type R2S1RPArray = Array U DIM4 Double -- (Z :. numOrientations :. numScales :. xLen :. yLen)
type R2S1T0Array = Array U DIM4 (Complex Double)  -- (Z :. (L.length theta0freqs) :. numOrientations :. xLen :. yLen )
type R2Z1T0Array = DFTArray DIM4 (Complex Double) -- (Z :. (L.length thetafreqs) :. (L.length theta0Freqs) :. xLen :. yLen )

type R2T0Array   = DFTArray DIM3 (Complex Double) -- (Z :. (L.length theta0freqs) :. xLen :. yLen)
type R2Z1Array   = DFTArray DIM3 (Complex Double) -- (Z :. (L.length thetafreqs) :. xLen :. yLen)

data R2S1RPPoint =
  R2S1RPPoint (Int, Int, Double, Double)
  deriving (Show, Read)

{-# INLINE repa2vec #-}
repa2vec :: (Shape dim, Unbox e, Storable e) => DFTArray dim e -> DFTArray dim e
repa2vec (RepaArray arr) =
  VectorArray (listOfShape . extent $ arr) . VU.convert . toUnboxed $ arr
  
{-# INLINE vec2repa #-}
vec2repa :: (Shape dim, Unbox e, Storable e) => DFTArray dim e -> DFTArray dim e
vec2repa (VectorArray xs vec) =
  RepaArray . fromUnboxed (shapeOfList xs) . VS.convert $ vec
