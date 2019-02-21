module FokkerPlanck.Types where

type ParticleIndex
   = ( Double -- x
     , Double -- y
     , Double -- theta
     , Double -- scale
     , Double -- theta_0
     , Double -- scale_0
      )

-- data ParticleIndex = R2S1
--   { getR2S1X     :: Double
--   , getR2S1Y     :: Double
--   , getR2S1Theta :: Double
--   } deriving (Show)
