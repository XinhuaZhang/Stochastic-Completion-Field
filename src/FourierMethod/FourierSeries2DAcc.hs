{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Strict              #-}
{-# LANGUAGE TypeOperators       #-}
module FourierMethod.FourierSeries2DAcc where

import           Data.Array.Accelerate              as A
import           Data.Array.Accelerate.Data.Complex as A

{-# INLINE harmonicAcc #-}
harmonicAcc ::
     (A.Floating e, Elt (Complex e), A.FromIntegral Int e)
  => Int
  -> e
  -> e
  -> e
  -> Acc (A.Vector (Int, Int))
  -> Acc (A.Array DIM3 (A.Complex e))
harmonicAcc numPoints period' delta' deltaFreq' vec =
  let period = constant period'
      delta = constant delta'
      deltaFreq = constant deltaFreq'
      center = constant $ div numPoints 2
      c = (-2) * A.pi / period * delta * deltaFreq
  in A.imap
       (\(unlift -> Z :. (_ :: Exp Int) :. x :. y) (unlift -> (xFreq, yFreq)) ->
          A.cis $
          c * A.fromIntegral (xFreq * (x - center) + yFreq * (y - center))) .
     A.replicate (constant (Z :. All :. numPoints :. numPoints)) $
     vec

{-# INLINE inverseHarmonicAcc #-}
inverseHarmonicAcc ::
     (A.Floating e, Elt (Complex e), A.FromIntegral Int e)
  => Int
  -> e
  -> e
  -> Acc (A.Vector (Int, Int))
  -> Acc (A.Array DIM3 (A.Complex e))
inverseHarmonicAcc numFreqs period' delta' vec =
  let period = constant period'
      delta = constant delta'
      center = constant $ div numFreqs 2
      c = 2 * A.pi / period * delta
  in A.imap
       (\(unlift -> Z :. xFreq :. yFreq :. (_ :: Exp Int)) (unlift -> (x, y)) ->
          A.cis $
          c * A.fromIntegral (x * (xFreq - center) + y * (yFreq - center))) .
     A.replicate (constant (Z :. numFreqs :. numFreqs :. All)) $
     vec
     
{-# INLINE inverseHarmonicAcc1 #-}
inverseHarmonicAcc1 ::
     (A.Floating e, Elt (Complex e), A.FromIntegral Int e)
  => Int
  -> e
  -> e
  -> Acc (A.Vector (Int, Int))
  -> Acc (A.Array DIM3 (A.Complex e))
inverseHarmonicAcc1 numFreqs period' delta' vec =
  let period = constant period'
      delta = constant delta'
      center = constant $ div numFreqs 2
      c = 2 * A.pi / period * delta
  in A.imap
       (\(unlift -> Z :. (_ :: Exp Int) :. xFreq :. yFreq) (unlift -> (x, y)) ->
          A.cis $
          c * A.fromIntegral (x * (xFreq - center) + y * (yFreq - center))) .
     A.replicate (constant (Z :. All :. numFreqs :. numFreqs)) $
     vec
