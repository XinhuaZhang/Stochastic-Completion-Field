{-# LANGUAGE BangPatterns  #-}
{-# LANGUAGE DeriveGeneric #-}
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
import           Control.Monad
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
import           Utils.Distribution

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
    Just d  -> genContVar d gen

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
  -> Double
  -> Particle
  -> DList Particle
  -> IO (DList Particle)
brownianMotion randomGen thetaDist scaleDist poissonDist deltaT maxRho maxR tao stdR2 particle xs = do
  deltaTheta <- generateRandomNumber randomGen thetaDist
  deltaScale <- generateRandomNumber randomGen scaleDist
  let newParticle@(Particle phi rho theta r v) = moveParticle deltaT particle
      ys =
        if r <= maxR && r > 1 / maxR && rho <= maxRho && rho > 1 / maxRho
          then DL.cons
                 (Particle phi rho theta r (1 - gaussian2DPolar rho stdR2))
                 xs
          else xs
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
           stdR2
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
  -> Double
  -> GenIO
  -> IO (DList Particle)
generatePath thetaDist scaleDist poissonDist maxRho maxR tao deltaT stdR2 randomGen =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    poissonDist
    deltaT
    maxRho
    maxR
    tao
    stdR2
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
  => GenIO ->  Double -> Maybe d -> Double -> Particle -> IO Particle
diffuseParticle' randomGen thetaSigma scaleDist deltaT (Particle phi rho theta r v) = do
  deltaScale <- generateRandomNumber randomGen scaleDist
  let r1 = r `scalePlus` deltaScale
  deltaTheta <-
    genContVar (normalDistr 0 (thetaSigma * sqrt (r1 * deltaT))) randomGen
  return $ Particle phi rho (theta `thetaPlus` deltaTheta) r1 v

brownianMotion' ::
     (Distribution d, ContGen d)
  => GenIO
  -> Double
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Particle
  -> DList Particle
  -> IO (DList Particle)
brownianMotion' randomGen thetaSigma scaleDist poissonLambda deltaT maxRho maxR tao stdR2 particle xs = do
  let newParticle@(Particle phi rho theta r v) = moveParticle deltaT particle
      ys =
        if r <= maxR && r > 1 / maxR && rho <= maxRho && rho >= 1 / maxRho
          then DL.cons (Particle phi rho theta r (1 - gaussian2DPolar rho stdR2)) xs
          else xs
  t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
  if t > exp (1 / (-tao / (r * deltaT)))
    then return ys
    else do
      diffusedParticle <-
        diffuseParticle' randomGen thetaSigma scaleDist deltaT newParticle
      brownianMotion'
        randomGen
        thetaSigma
        scaleDist
        poissonLambda
        deltaT
        maxRho
        maxR
        tao
        stdR2
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
  -> Double
  -> Double
  -> GenIO
  -> IO (DList Particle)
generatePath' thetaSigma scaleDist poissonLambda maxRho maxR tao deltaT stdR2 randomGen = do
  brownianMotion'
    randomGen
    thetaSigma
    scaleDist
    poissonLambda
    deltaT
    maxRho
    maxR
    tao
    stdR2
    (Particle 0 1 0 1 1)
    DL.empty
