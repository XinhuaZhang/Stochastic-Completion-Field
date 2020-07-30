{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
module FourierPinwheel.Array where

data FPArray vector = FPArray
  { getFPArrayNumXFreq     :: Int
  , getFPArrayNumYFreq     :: Int
  , getFPArrayNumRFreq     :: Int
  , getFPArrayNumThetaFreq :: Int
  , getFPArrayNumRhoFreq   :: Int
  , getFPArrayNumPhiFreq   :: Int
  , getFPArray             :: [vector]
  }

