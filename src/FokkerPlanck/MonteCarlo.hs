{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.MonteCarlo
  ( module FokkerPlanck.Types
  , solveMonteCarloR2S1
  , solveMonteCarloR2S1RP
  ) where

import           Array.UnboxedArray              as UA
import           Control.Arrow
import           Control.Monad                   as M
import           Control.Monad.Parallel          as MP
import           Data.Array                      as Arr
import           Data.Array.Repa                 as R
import           Data.Complex
import           Data.DList                      as DL
import           Data.List                       as L
import           Data.Vector.Unboxed             as VU
import           FokkerPlanck.Types
import           Statistics.Distribution
import           Statistics.Distribution.Normal
import           Statistics.Distribution.Uniform
import           System.Random
import           System.Random.MWC
import           Utils.Coordinates
import           Utils.Parallel

{-# INLINE thetaPlus #-}
thetaPlus :: Double -> Double -> Double
thetaPlus x y
  | z < 0 =  z + a
  | z >= a = z - a
  | otherwise = z
  where
    !z = x + y
    a = 2 * pi

{-# INLINE thetaCheck #-}
thetaCheck :: Double -> Double
thetaCheck !theta =
  if theta < 0
    then theta + 2 * pi
    else if theta >= 2 * pi
           then theta - 2 * pi
           else theta

{-# INLINE scalePlus #-}
scalePlus :: Double -> Double -> Double -> (Double, Bool)
scalePlus maxScale x y
  | z < 0 = (-z, True)
  | z >= maxScale = (maxScale, False)
  | otherwise = (z, False)
  where
    !z = x + y

{-# INLINE generatePath #-}
generatePath ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePath randomGen thetaDist' scaleDist' maxScale tao numSteps init@(x', y', theta', scale', theta0', scale0') =
  go
    numSteps
    (x', y', thetaCheck theta', scale', thetaCheck theta0', scale0')
    DL.empty
  where
    go 0 _ xs = return xs
    go n z@(!x, !y, theta, scale, theta0, scale0) xs = do
      deltaTheta <-
        case thetaDist' of
          Nothing -> return 0
          Just thetaDist -> genContVar thetaDist randomGen
      deltaScale <-
        case scaleDist' of
          Nothing -> return 0
          Just scaleDist -> genContVar scaleDist randomGen
      let (newScale, flag) = scalePlus maxScale scale deltaScale
          newTheta =
            if flag
              then (theta `thetaPlus` (deltaTheta + pi))
              else (theta `thetaPlus` deltaTheta)
          !newX = x + newScale * cos theta
          !newY = y + newScale * sin theta
          newIndex = (newX, newY, newTheta, newScale, theta0, scale0)
      t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
      if t < (1 - exp ((-1) / tao))
        then go 0 newIndex (DL.cons newIndex xs)
        else go (n - 1) newIndex (DL.cons newIndex xs)


{-# INLINE countR2S1 #-}
countR2S1 :: Int -> Int -> Int -> [DList ParticleIndex] -> (Int, VU.Vector Int)
countR2S1 xLen yLen numOrientations xs =
  let !xShift = div xLen 2
      (!xMin, !xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      (!yMin, !yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        VU.fromList .
        L.filter
          (\(x, y, _) ->
             (x >= xMin) && (x <= xMax) && (y >= yMin) && (y <= yMax)) .
        DL.toList .
        DL.concat .
        L.map (DL.map (\(x, y, t, s, t0, s0) -> (round x, round y, t))) $
        xs
      numTrajectories = VU.length ys
      arr =
        UA.accumulate
          (+)
          (0 :: Int)
          ((xMin, yMin, 0), (xMax, yMax, numOrientations - 1)) .
        VU.map (\(x, y, t) -> ((x, y, floor (t / deltaTheta)), 1)) $
        ys
   in (numTrajectories, toUnboxedVector arr)

{-# INLINE solveMonteCarloR2S1 #-}
solveMonteCarloR2S1 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (R.Array U DIM3 Double)
solveMonteCarloR2S1 numGen numTrails xLen yLen numOrientations thetaSigma tao numSteps (_, _, _, _, t0, _) = do
  gens <- M.replicateM numGen createSystemRandom
  let thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 0
  xs <-
    MP.mapM
      (\gen ->
         M.replicateM
           (div numTrails numGen)
           (generatePath
              gen
              thetaDist
              scaleDist
              10
              tao
              numSteps
              (0, 0, t0, 1, t0, 1)))
      gens
  let (ys, zs) =
        L.unzip $ parMap rdeepseq (countR2S1 xLen yLen numOrientations) xs
      totalNum = fromIntegral $ L.sum ys :: Double
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  print totalNum
  print $ L.length zs
  return .
    fromUnboxed (Z :. xLen :. yLen :. numOrientations) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec

{-# INLINE countR2S1RP #-}
countR2S1RP ::
     Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> [DList ParticleIndex]
  -> (Int, VU.Vector Int)
countR2S1RP xLen yLen numOrientations numScales maxScale xs =
  let !xShift = div xLen 2
      (!xMin, !xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      (!yMin, !yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      !deltaScale = maxScale / (fromIntegral numScales)
      ys =
        VU.fromList .
        L.filter
          (\(x, y, _, s) ->
             (x >= xMin) &&
             (x <= xMax) && (y >= yMin) && (y <= yMax) -- && (s /= 0)
          ) .
        DL.toList .
        DL.concat .
        L.map (DL.map (\(x, y, t, s, t0, s0) -> (round x, round y, t, s))) $
        xs
      numTrajectories = VU.length ys
      arr =
        UA.accumulate
          (+)
          0
          ((xMin, yMin, 0, 0), (xMax, yMax, numOrientations - 1, numScales - 1)) .
        (VU.map
           (\(x, y, t, s) ->
              ( ( x
                , y
                , floor $ t / deltaTheta
                , (let s1 = floor $ s / deltaScale
                    in if s1 == numScales
                         then numScales - 1
                         else s1))
              , 1))) $
        ys
   in (numTrajectories, toUnboxedVector arr)

-- -- total number of trails = numGen * numTrails
{-# INLINE solveMonteCarloR2S1RP #-}
solveMonteCarloR2S1RP ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> ParticleIndex
  -> IO (R.Array U DIM4 Double)
solveMonteCarloR2S1RP numGen numTrails xLen yLen numOrientations numScales thetaSigma scaleSigma maxScale tao numSteps (_, _, _, _, t0, s0) = do
  gens <- M.replicateM numGen createSystemRandom
  let thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 scaleSigma
  xs <-
    MP.mapM
      (\gen ->
         M.replicateM
           (div numTrails numGen)
           (generatePath
              gen
              thetaDist
              scaleDist
              maxScale
              tao
              numSteps
              (0, 0, t0, s0, t0, s0)))
      gens
  let (ys, zs) =
        L.unzip $
        parMap
          rdeepseq
          (countR2S1RP xLen yLen numOrientations numScales maxScale)
          xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  return .
    fromUnboxed (Z :. xLen :. yLen :. numOrientations :. numScales) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec
