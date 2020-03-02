{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.BrownianMotion
  ( Particle(..)
  , generatePath
  , moveParticle
  ) where

import           Control.DeepSeq
import           Data.DList                          as DL
import           Data.Ix
import           Data.List                           as L
import           GHC.Generics                        (Generic)
import           Statistics.Distribution
import           Statistics.Distribution.Exponential
import           Statistics.Distribution.Normal
import           Statistics.Distribution.Uniform
import           System.Random.MWC
import           Text.Printf

data Particle =
  Particle {-# UNPACK #-}!Double -- \phi
           {-# UNPACK #-}!Double -- \rho
           {-# UNPACK #-}!Double -- \theta
           {-# UNPACK #-}!Double -- r
  deriving (Show,Generic)

instance NFData Particle

{-# INLINE getParticleRho #-}
getParticleRho :: Particle -> Double
getParticleRho (Particle _ rho _ _ ) = rho

{-# INLINE thetaCheck #-}
thetaCheck :: Double -> Double
thetaCheck theta =
  if theta >= -pi && theta < pi
    then theta
    else error (printf "thetaCheck: %.2f is out of boundary (0,2pi)\n" theta)

{-# INLINE scaleCheck #-}
scaleCheck :: Double -> Double -> Double
scaleCheck maxScale scale =
  if scale >= (1 / maxScale) && scale <= maxScale
    then scale
    else error
           (printf
              "scaleCheck: %.2f is out of boundary (0,%.2f)\n"
              scale
              maxScale)

{-# INLINE thetaPlus #-}
thetaPlus :: Double -> Double -> Double
thetaPlus !x !y
  | z <= -pi = z + a
  | z > pi = z - a
  | otherwise = z
  where
    !z = x + y
    !a = 2 * pi
         
{-# INLINE scalePlus #-}
scalePlus :: Double -> Double -> Double  
scalePlus x delta = x * exp delta 

{-# INLINE scalePlusPeriodic #-}
scalePlusPeriodic :: Double -> Double -> Double -> Double
scalePlusPeriodic !logMaxScale !delta !x
  | z >= m = (exp (z - 2 * m))
  | z < -m = (exp (2 * m + z))
  | otherwise = (exp z)
  where
    !m = logMaxScale
    !z = (log x) + delta

-- {-# INLINE scalePlus #-}
-- scalePlus :: Double -> Double -> Double -> Double
-- scalePlus !maxScale !x !y
--   | z < -maxScale = z + 2 * maxScale
--   | z >= maxScale = z - 2 * maxScale
--   | otherwise = z
--   where
--     !z = x + y

{-# INLINE rhoCutoff #-}
rhoCutoff :: Double -> Double
rhoCutoff !rho =
  if rho < 0
    then 0
    else rho

-- {-# INLINE scalePlus' #-}
-- scalePlus' :: Double -> Double -> Double
-- scalePlus' !maxScale !z
--   | z < -maxScale = z + 2 * maxScale
--   | z >= maxScale = z - 2 * maxScale
--   | otherwise = z

{-# INLINE generateRandomNumber #-}
generateRandomNumber ::
     (Distribution d, ContGen d) => GenIO -> Maybe d -> IO Double
generateRandomNumber !gen dist =
  case dist of
    Nothing -> return 0
    Just d -> genContVar d gen

{-# INLINE moveParticle #-}
moveParticle :: Particle -> Particle
moveParticle (Particle phi rho theta r) =
  let !x = theta - phi
      !cosX = cos x
      !newPhi = phi `thetaPlus` atan2 (r * sin x) (rho + r * cosX)
      !newRho = sqrt $ (rho * rho + r * r + 2 * r * rho * cosX)
  in Particle newPhi newRho theta r
  
-- {-# INLINE moveParticle #-}
-- moveParticle :: Double -> Particle -> Particle
-- moveParticle !maxScale (Particle phi rho theta r) =
--   let !newPhi = phi `thetaPlus` theta
--       !newRho = scalePlus maxScale rho r
--   in Particle newPhi newRho theta r

{-# INLINE diffuseParticle #-}
diffuseParticle :: Double -> Double -> Double -> Particle -> Particle
diffuseParticle !deltaTheta !deltaScale !maxScale (Particle phi rho theta r) =
  Particle
    phi
    rho
    (theta `thetaPlus` deltaTheta)
    (r `scalePlus` deltaScale)
    -- (scalePlusPeriodic (log maxScale) deltaScale r)

brownianMotion ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Particle
  -> DList Particle
  -> IO (DList Particle)
brownianMotion !randomGen !thetaDist !scaleDist !maxScale !tao !particle !xs = do
  deltaTheta <- generateRandomNumber randomGen thetaDist
  deltaScale <- generateRandomNumber randomGen scaleDist
  let !newParticle@(Particle !phi !rho !theta !r) =
            diffuseParticle deltaTheta deltaScale maxScale . moveParticle  $  particle
      ys = -- DL.cons (Particle phi (rho) theta (r)) xs
        if rho > 0 -- && r > 0
          then DL.cons (Particle phi (rho) theta (r)) xs
          else xs
  t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
  if t > exp (1 / (-tao))
    then return ys
    else brownianMotion
           randomGen
           thetaDist
           scaleDist
           maxScale
           tao
           (newParticle)
           ys

{-# INLINE generatePath #-}
generatePath ::
     (Distribution d, ContGen d)
  => Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> GenIO
  -> IO (DList Particle)
generatePath thetaDist scaleDist maxScale tao initScale randomGen = do
  -- deltaTheta <- generateRandomNumber randomGen thetaDist
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    maxScale
    tao
    (Particle 0 1 0 initScale)
    DL.empty
    -- (DL.fromList [(Particle 0 0 0 1)])

