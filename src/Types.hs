{-# LANGUAGE TypeOperators #-}
module Types where

import           Data.Array.Repa
import           Data.Complex
                
type DIM6 = DIM5 :. Int

type R2S1Array   = Array U DIM3 (Complex Double) -- (Z :. numOrientations :. xLen :. yLen)
type R2S1RPArray = Array U DIM4 Double            -- (Z :. numOrientations :. numScales :. xLen :. yLen)
type R2S1T0Array = Array U DIM4 (Complex Double)  -- (Z :. numOrientations :. (L.length theta0freqs)  :. xLen :. yLen )
type R2Z1T0Array = Array U DIM4 (Complex Double) -- (Z :. (L.length thetafreqs) :. (L.length theta0Freqs) :. xLen :. yLen )


type R2S1RPT0S0Array = Array U DIM6 (Complex Double) 
-- (Z :. numOrientations :. numScales :. (L.length theta0freqs) :. (L.length scale0freqs)  :. xLen :. yLen )
type R2Z2T0S0Array = Array U DIM6 (Complex Double) 
-- (Z :. (L.length thetafreqs)  :. (L.length scalefreqs) :. (L.length theta0Freqs) :. (L.length scale0Freqs) :. xLen :. yLen )

type R2T0Array   = Array U DIM3 (Complex Double) -- (Z :. (L.length theta0freqs) :. xLen :. yLen)
type R2Z1Array   = Array U DIM3 (Complex Double) -- (Z :. (L.length thetafreqs) :. xLen :. yLen)

type R2T0S0Array   = Array U DIM4 (Complex Double) -- (Z :. (L.length theta0freqs) :. (L.length scale0freqs)  :. xLen :. yLen)

data R2S1RPPoint =
  R2S1RPPoint (Int, Int, Double, Double)
  deriving (Show, Read)
