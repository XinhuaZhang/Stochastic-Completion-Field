module FourierPinwheel.Hypergeo1F1 where

import           Data.Complex
import           Data.List              as L
import           Math.Sequence.Converge

{-# INLINE hypergeom #-}
hypergeom ::
     (Eq a, RealFloat a, Enum a)
  => Complex a
  -> Complex a
  -> Complex a
  -> Complex a
hypergeom a b z =
  converge .
  scanl (+) 0 .
  scanl (*) 1 .
  L.map (\i -> (a + ((i - 1) :+ 0)) / (b + ((i - 1) :+ 0)) * z / (i :+ 0)) $
  [1 ..]
