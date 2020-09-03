{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.BrownianMotion
  ( Particle(..)
  , generatePath
  , moveParticle
  , thetaPlus
  , scalePlus
  , generateRandomNumber
  , generatePath'
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
import Control.Monad

data Particle =
  Particle {-# UNPACK #-}!Double -- \phi
           {-# UNPACK #-}!Double -- \rho
           {-# UNPACK #-}!Double -- \theta
           {-# UNPACK #-}!Double -- r
           {-# UNPACK #-}!Double -- weight
  deriving (Show,Generic)

instance NFData Particle

{-# INLINE getParticleRho #-}
getParticleRho :: Particle -> Double
getParticleRho (Particle _ rho _ _ _ ) = rho


thetaCheck :: Double -> Double
thetaCheck theta =
  if theta < -pi
    then thetaCheck (theta + 2 * pi)
    else if theta >= pi
           then thetaCheck (theta - 2 * pi)
           else theta

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
thetaPlus !x !y =
  let z = x + y
   in thetaCheck z
         
{-# INLINE scalePlus #-}
scalePlus :: Double -> Double -> Double  
scalePlus x delta = x * exp delta 
          
{-# INLINE scalePlusPeriodic #-}
scalePlusPeriodic :: Double -> Double -> Double -> Double
scalePlusPeriodic !logMaxScale !delta !x
  | z >= m = exp (z - 2 * m)
  | z < -m = exp (2 * m + z)
  | otherwise = exp z
  where
    !m = logMaxScale
    !z = log x + delta
         
{-# INLINE rhoCutoff #-}
rhoCutoff :: Double -> Double
rhoCutoff !rho =
  if rho < 0
    then 0
    else rho

{-# INLINE generateRandomNumber #-}
generateRandomNumber ::
     (Distribution d, ContGen d) => GenIO -> Maybe d -> IO Double
generateRandomNumber !gen dist =
  case dist of
    Nothing -> return 0
    Just d -> genContVar d gen

{-# INLINE moveParticle #-}
moveParticle :: Double -> Particle -> Particle
moveParticle deltaT (Particle phi rho theta r0 v) =
  let !r = r0 * deltaT
      !x = theta - phi
      !cosX = cos x
      !newPhi = phi `thetaPlus` atan2 (r * sin x) (rho + r * cosX)
      !newRho = sqrt (rho * rho + r * r + 2 * r * rho * cosX)
  in Particle newPhi newRho theta r0 v

{-# INLINE diffuseParticle #-}
diffuseParticle :: Double -> Double -> Double -> Particle -> Particle
diffuseParticle !deltaTheta !deltaScale maxScale (Particle phi rho theta r v) =
  Particle
    phi
    rho
    (theta `thetaPlus` deltaTheta)
    (r `scalePlus` deltaScale)
    v

brownianMotion ::
     (Distribution d1, ContGen d1, Distribution d2, ContGen d2, Distribution d3)
  => GenIO
  -> Maybe d1
  -> Maybe d2
  -> Maybe d3
  -> Double
  -> Double
  -> Double
  -> Double
  -> Particle
  -> DList Particle
  -> IO (DList Particle)
brownianMotion randomGen thetaDist scaleDist poissonDist deltaT maxRho maxR tao particle xs = do
  deltaTheta <- generateRandomNumber randomGen thetaDist
  deltaScale <- generateRandomNumber randomGen scaleDist
  let newParticle@(Particle phi rho theta r v) = moveParticle deltaT particle
      -- phi1 = phi `thetaPlus` (-theta)
      -- rho1 = rho / r
      -- theta1 = -theta
      -- r1 = 1 / r
      ys =
        if r <= maxR && r > 1 / maxR   && rho <= maxRho && rho >= 1
          then DL.cons newParticle xs
          else xs
      -- ys =
      --   if rho1 <= maxR && rho1 > 1
      --     then DL.cons (Particle phi1 rho1 theta1 r1 v) ys'
      --     else ys'
  t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
  if t > exp (1 / (-tao))
    then return ys
    else brownianMotion
           randomGen
           thetaDist
           scaleDist
           poissonDist
           deltaT
           maxRho
           maxR
           tao
           (diffuseParticle deltaTheta deltaScale maxR newParticle)
           ys
           
{-# INLINE generatePath #-}
generatePath ::
     (Distribution d1, ContGen d1, Distribution d2, ContGen d2, Distribution d3)
  => Maybe d1
  -> Maybe d2
  -> Maybe d3
  -> Double
  -> Double
  -> Double
  -> Double
  -> GenIO
  -> IO (DList Particle)
generatePath thetaDist scaleDist poissonDist maxRho maxR tao deltaT randomGen =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    poissonDist
    deltaT
    maxRho
    maxR
    tao
    (Particle 0 1 0 1 1)
    DL.empty

{-# INLINE moveParticle' #-}
moveParticle' ::  Particle -> Particle
moveParticle' (Particle phi rho theta r v) =
  let !x = theta - phi
      !cosX = cos x
      !newPhi = phi `thetaPlus` atan2 (r * sin x) (rho + r * cosX)
      !newRho = sqrt (rho * rho + r * r + 2 * r * rho * cosX)
   in Particle newPhi newRho theta r v

{-# INLINE diffuseParticle' #-}
diffuseParticle' :: (Distribution d, ContGen d)
  => GenIO ->  Double -> Maybe d -> Particle -> IO Particle
diffuseParticle' randomGen thetaSigma scaleDist (Particle phi rho theta r v) = do
  deltaScale <- generateRandomNumber randomGen scaleDist
  let r1 = r `scalePlus` deltaScale
      thetaDist = normalDistrE 0 (thetaSigma * sqrt r1)
  deltaTheta <- generateRandomNumber randomGen thetaDist
  return $
    Particle phi rho (theta `thetaPlus` deltaTheta) r1 v

brownianMotion' ::
     (Distribution d, ContGen d)
  => GenIO
  -> Double
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> Particle
  -> DList Particle
  -> IO (DList Particle)
brownianMotion' randomGen thetaSigma scaleDist poissonLambda maxScale tao particle@(Particle _ _ _ r0 _) xs = do
  let newParticle@(Particle phi rho theta r v) = moveParticle' particle
      ys =
        if r <= maxScale && r > 1 / maxScale && rho <= maxScale && rho > 1
          then DL.cons newParticle xs
          else xs
  t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
  if t > exp (1 / (-tao / r0))
    then return ys
    else do
      diffusedParticle <-
        diffuseParticle' randomGen thetaSigma scaleDist newParticle
      brownianMotion'
        randomGen
        thetaSigma
        scaleDist
        poissonLambda
        maxScale
        tao
        diffusedParticle
        ys


{-# INLINE generatePath' #-}
generatePath' ::
     (Distribution d, ContGen d)
  => Double
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> Double
  -> GenIO
  -> IO (DList Particle)
generatePath' thetaSigma scaleDist poissonLambda tao deltaT maxScale randomGen = do
  brownianMotion'
    randomGen
    thetaSigma
    scaleDist
    poissonLambda
    maxScale
    tao
    (Particle 0 0 0 deltaT 1)
    DL.empty
