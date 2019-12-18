{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.BrownianMotion 
  ( generatePath
  , generatePathRaidalRounded
  , generatePathRaidalAntiAliasing
  ) where

import           Data.DList                          as DL
import           Data.Ix
import           Data.List                           as L
import           FokkerPlanck.Types
import           Statistics.Distribution
import           Statistics.Distribution.Exponential
import           Statistics.Distribution.Normal
import           Statistics.Distribution.Uniform
import           System.Random.MWC
import           Text.Printf

{-# INLINE thetaCheck #-}
thetaCheck :: Double -> Double
thetaCheck theta =
  if theta >= 0 && theta <= (2 * pi)
    then theta
    else error (printf "thetaCheck: %.2f is out of boundary (0,2pi)\n" theta)
    
{-# INLINE scaleCheck #-}
scaleCheck :: Double -> Double -> Double
scaleCheck maxScale scale =
  if scale >= 0 && scale <= maxScale
    then scale
    else error
           (printf
              "scaleCheck: %.2f is out of boundary (0,%.2f)\n"
              scale
              maxScale)

{-# INLINE thetaPlus #-}
thetaPlus :: Double -> Double -> Double
thetaPlus x y
  | z < 0 = z + a
  | z >= a = z - a
  | otherwise = z
  where
    z = x + y
    a = 2 * pi

{-# INLINE scalePlus #-}
scalePlus :: Double -> Double -> Double -> Double
scalePlus maxScale x y
  | z < 0 = 0
  | z >= maxScale = maxScale
  | otherwise = z
  where
    z = x + y

{-# INLINE scalePlusPeriodic #-}
scalePlusPeriodic :: Double -> Double -> Double -> Double
scalePlusPeriodic maxScale x y
  | z < 0 = maxScale + z
  | z >= maxScale = z - maxScale
  | otherwise = z
  where
    z = x + y

{-# INLINE generateRandomNumber #-}
generateRandomNumber ::
     (Distribution d, ContGen d) => GenIO -> Maybe d -> IO Double
generateRandomNumber gen dist =
  case dist of
    Nothing -> return 0
    Just d -> genContVar d gen

{-# INLINE brownianMotion #-}
brownianMotion ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> (Double -> Double -> Bool)
  -> ParticleIndex
  -> DList ParticleIndex
  -> IO (DList ParticleIndex)
brownianMotion randomGen thetaDist scaleDist maxScale tao deltaDist checkRange (!x, !y, !theta, !scale, theta0, scale0) xs = do
  deltaTheta <- generateRandomNumber randomGen thetaDist
  deltaScale <- generateRandomNumber randomGen scaleDist
  let newTheta = theta `thetaPlus` deltaTheta
      newScale = scalePlus maxScale scale deltaScale
      cosTheta = cos theta
      sinTheta = sin theta
      newX = x + scale * cosTheta
      newY = y + scale * sinTheta
      newIndex = (newX, newY, newTheta, newScale, theta0, scale0)
      ys =
        if checkRange newX newY
          then DL.cons newIndex xs
          else xs
      iDeltaDist i = deltaDist * fromIntegral i
      zs =
        DL.append ys .
        DL.fromList . L.filter (\(a, b, _, _, _, _) -> checkRange a b) $
        [ ( x + iDeltaDist i * cosTheta
          , y + iDeltaDist i * sinTheta
          , theta
          , iDeltaDist i
          , theta0
          , scale0)
        | i <- [1 .. (floor $ scale / deltaDist) - 1]
        ]
  t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
  if t < (1 - exp ((-1) / tao))
    then return zs
    else brownianMotion
           randomGen
           thetaDist
           scaleDist
           maxScale
           tao
           deltaDist
           checkRange
           newIndex
           zs

{-# INLINE generatePath #-}
generatePath ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> (Int, Int)
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePath randomGen thetaDist scaleDist maxScale tao deltaDist xRange yRange init@(x', y', theta', scale', theta0', scale0') =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    maxScale
    tao
    deltaDist
    (\a b ->
       (inRange xRange (round a)) &&
       (inRange yRange (round b)) && (round a /= 0 || round b /= 0))
    ( x'
    , y'
    , thetaCheck theta'
    , scaleCheck maxScale scale'
    , thetaCheck theta0'
    , scaleCheck maxScale scale0')
    DL.empty

-- Only simulate the radial part
{-# INLINE generatePathRaidalRounded #-}
generatePathRaidalRounded ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePathRaidalRounded randomGen thetaDist scaleDist maxScale tao deltaDist rRange init@(x', y', theta', scale', theta0', scale0') =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    maxScale
    tao
    deltaDist
    (\a b -> (inRange rRange (round a)) && (round b == 0) && (round a /= 0))
    ( x'
    , y'
    , thetaCheck theta'
    , scaleCheck maxScale scale'
    , thetaCheck theta0'
    , scaleCheck maxScale scale0')
    DL.empty
    

{-# INLINE generatePathRaidalAntiAliasing #-}
generatePathRaidalAntiAliasing ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePathRaidalAntiAliasing randomGen thetaDist scaleDist maxScale tao deltaDist rRange init@(x', y', theta', scale', theta0', scale0') =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    maxScale
    tao
    deltaDist
    (\x y ->
       x > fromIntegral (rMin - 1) &&
       x < fromIntegral (rMax + 1) && y < 1 && y > -1)
    ( x'
    , y'
    , thetaCheck theta'
    , scaleCheck maxScale scale'
    , thetaCheck theta0'
    , scaleCheck maxScale scale0')
    DL.empty
  where
    (rMin, rMax) = rRange


{-# INLINE generatePathRaidalGaussian #-}
generatePathRaidalGaussian ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> (Int, Int)
  -> Double
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePathRaidalGaussian randomGen thetaDist scaleDist maxScale tao deltaDist rRange gaussianRadius init@(x', y', theta', scale', theta0', scale0') =
  brownianMotion
    randomGen
    thetaDist
    scaleDist
    maxScale
    tao
    deltaDist
    (\x y ->
       x >= -gaussianRadius &&
       x <= gaussianRadius + fromIntegral rMax &&
       y <= gaussianRadius && y >= -gaussianRadius)
    ( x'
    , y'
    , thetaCheck theta'
    , scaleCheck maxScale scale'
    , thetaCheck theta0'
    , scaleCheck maxScale scale0')
    DL.empty
  where
    (rMin, rMax) = rRange
