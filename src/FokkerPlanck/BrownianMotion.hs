{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.BrownianMotion
  ( Particle(..)
  , generatePath
  , moveParticle
  ) where

import           Data.DList                          as DL
import           Data.Ix
import           Data.List                           as L
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
  deriving (Show)

{-# INLINE getParticleRho #-}
getParticleRho :: Particle -> Double
getParticleRho (Particle _ rho _ _ ) = rho

{-# INLINE thetaCheck #-}
thetaCheck :: Double -> Double
thetaCheck theta =
  if theta >= 0 && theta <= (2 * pi)
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
  | z < 0 = z + a
  | z >= a = z - a
  | otherwise = z
  where
    !z = x + y
    !a = 2 * pi

{-# INLINE scalePlusPeriodic #-}
scalePlusPeriodic :: Double -> Double -> Double -> Double
scalePlusPeriodic !maxScale !x !y
  | z < m = maxScale - (m - z)
  | z >= maxScale = z - maxScale + m
  | otherwise = z
  where
    !m = 1 / maxScale
    !z = x + y

{-# INLINE scalePlusPeriodic' #-}
scalePlusPeriodic' :: Double -> Double -> Double  
scalePlusPeriodic' !maxScale !z
  | z < m = maxScale - (m - z)
  | z >= maxScale = z - maxScale + m
  | otherwise = z
  where
    m = 1 / maxScale

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
  let x = theta - phi
      cosX = cos x
      newPhi = phi `thetaPlus` atan2 (r * sin x) (rho + r * cosX)
      newRho = sqrt (rho * rho + r * r + 2 * r * rho * cosX)
  in Particle newPhi newRho theta r

{-# INLINE diffuseParticle #-}
diffuseParticle :: Double -> Double -> Double -> Particle -> Particle
diffuseParticle !deltaTheta !deltaScale !maxScale (Particle phi rho theta r) =
  Particle
    phi
    rho
    (theta `thetaPlus` deltaTheta)
    (scalePlusPeriodic maxScale r deltaScale)

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
  let !newParticle@(Particle _ !rho _ _) =
        diffuseParticle deltaTheta deltaScale maxScale . moveParticle $ particle
      ys =
        if 1 / maxScale <= rho && rho <= maxScale
          then DL.cons newParticle xs
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
           newParticle
           ys

{-# INLINE generatePath #-}
generatePath ::
     (Distribution d, ContGen d)
  => Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> GenIO
  -> IO (DList Particle)
generatePath thetaDist scaleDist maxScale tao randomGen =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    maxScale
    tao
    (Particle 0 0 0 1)
    DL.empty
