{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.MonteCarlo
  ( module FokkerPlanck.Types
  , solveMonteCarloR2S1
  , solveMonteCarloR2S1RP
  , solveMonteCarloR2S1T0
  , solveMonteCarloR2S1T0T
  ) where

import           Array.UnboxedArray              as UA
import           Control.Arrow
import           Control.Monad                   as M
import           Control.Monad.Parallel          as MP
import           Data.Array                      as Arr
import           Data.Array.Repa                 as R
import           Data.Complex
import           Data.DList                      as DL
import           Data.Ix
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
  -> (Int, Int)
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePath randomGen thetaDist' scaleDist' maxScale tao numSteps xRange yRange init@(x', y', theta', scale', theta0', scale0') =
  go
    numSteps
    (x', y', thetaCheck theta', scale', thetaCheck theta0', scale0')
    DL.empty
  where
    go 0 _ xs = return xs
    go n z@(!x, !y, !theta, !scale, theta0, scale0) xs = do
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
          !newX = x + scale * cos theta
          !newY = y + scale * sin theta
          newIndex = (newX, newY, newTheta, newScale, theta0, scale0)
          ys =
            if (inRange xRange (round newX)) && (inRange yRange (round newY))
              then DL.cons newIndex xs
              else xs
      t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
      if t < (1 - exp ((-1) / tao))
        then go 0 newIndex ys
        else go (n - 1) newIndex ys

-- R2S1

{-# INLINE countR2S1 #-}
countR2S1 ::
     (Int, Int)
  -> (Int, Int)
  -> Int
  -> [DList ParticleIndex]
  -> (Int, VU.Vector Int)
countR2S1 (xMin, xMax) (yMin, yMax) numOrientations xs =
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        L.map (\(t, x, y) -> ((floor $ t / deltaTheta, round x, round y), 1)) .
        DL.toList .
        DL.concat . L.map (DL.map (\(x, y, t, s, t0, s0) -> ((t, x, y)))) $
        xs
      numTrajectories = L.length ys
      arr =
        UA.accum (+) 0 ((0, xMin, yMin), (numOrientations - 1, xMax, yMax)) ys
   in (numTrajectories, toUnboxedVector arr)

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
  let !xShift = div xLen 2
      xRange =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      yRange =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      thetaDist = normalDistrE 0 thetaSigma
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
              xRange
              yRange
              (0, 0, t0, 1, t0, 1)))
      gens
  let (ys, zs) =
        L.unzip $ parMap rdeepseq (countR2S1 xRange yRange numOrientations) xs
      totalNum = fromIntegral $ L.sum ys :: Double
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  print totalNum
  print $ L.length zs
  return .
    fromUnboxed (Z :. numOrientations :. xLen :. yLen) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec

-- R2S1RP

{-# INLINE countR2S1RP #-}
countR2S1RP ::
     (Int, Int)
  -> (Int, Int)
  -> Int
  -> Int
  -> Double
  -> [DList ParticleIndex]
  -> (Int, VU.Vector Int)
countR2S1RP (xMin, xMax) (yMin, yMax) numOrientations numScales maxScale xs =
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      !deltaScale = maxScale / (fromIntegral numScales)
      ys =
        L.map (\(x, y, t, s, t0, s0) -> (t, s, round x, round y)) .
        DL.toList . DL.concat $
        xs
      numTrajectories = L.length ys
      arr =
        UA.accum
          (+)
          0
          ((0, 0, xMin, yMin), (numOrientations - 1, numScales - 1, xMax, yMax)) .
        (L.map
           (\(t, s, x, y) ->
              ( ( floor $ t / deltaTheta
                , (let s1 = floor $ s / deltaScale
                    in if s1 == numScales
                         then numScales - 1
                         else s1)
                , x
                , y)
              , 1))) $
        ys
   in (numTrajectories, toUnboxedVector arr)

-- -- total number of trails = numGen * numTrails
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
  let !xShift = div xLen 2
      xRange =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      yRange =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      thetaDist = normalDistrE 0 thetaSigma
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
              xRange
              yRange
              (0, 0, t0, s0, t0, s0)))
      gens
  let (ys, zs) =
        L.unzip $
        parMap
          rdeepseq
          (countR2S1RP xRange yRange numOrientations numScales maxScale)
          xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  return .
    fromUnboxed (Z :. numOrientations :. numScales :. xLen :. yLen) .
    VU.map (\x -> fromIntegral x / totalNum) $
    totalNumVec

-- R2S1 \theta_0 is represented in the frequency domain
{-# INLINE countR2S1T0 #-}
countR2S1T0 ::
     (Int, Int)
  -> (Int, Int)
  -> Int
  -> [Double]
  -> [DList ParticleIndex]
  -> (Int, VU.Vector (Complex Double))
countR2S1T0 (xMin, xMax) (yMin, yMax) numOrientations freqs xs =
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        DL.toList . DL.concat $
        L.zipWith
          (\freq i ->
             DL.map
               (\(x, y, t, _, t0, _) ->
                  ( (floor (t / deltaTheta), round x, round y, i)
                  , exp (0 :+ (freq * t0)))) .
             DL.concat $
             xs)
          freqs
          [1 ..]
      numTrajectories = L.length ys
      arr =
        UA.accum
          (+)
          0
          ( (0, xMin, yMin, 1)
          , (numOrientations - 1, xMax, yMax, L.length freqs))
          ys
   in (numTrajectories, toUnboxedVector arr)

-- (Z :. numOrientations :. xLen :. yLen :. (L.length freqs))
solveMonteCarloR2S1T0 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> [Double]
  -> ParticleIndex
  -> IO (R.Array U DIM4 (Complex Double))
solveMonteCarloR2S1T0 numGen numTrails xLen yLen numOrientations thetaSigma tao numSteps theta0freqs _ = do
  gens <- M.replicateM numGen createSystemRandom
  let !xShift = div xLen 2
      xRange =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      yRange =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 0
  xs <-
    MP.mapM
      (\gen ->
         M.replicateM
           (div numTrails numGen)
           (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
               generatePath
                 gen
                 thetaDist
                 scaleDist
                 10
                 tao
                 numSteps
                 xRange
                 yRange
                 (0, 0, t0, 1, t0, 1)))
      gens
  let (ys, zs) =
        L.unzip $
        parMap
          rdeepseq
          (countR2S1T0 xRange yRange numOrientations theta0freqs)
          xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  print totalNum
  print $ L.length zs
  return .
    fromUnboxed (Z :. numOrientations :. xLen :. yLen :. (L.length theta0freqs)) .
    VU.map (\x -> x / totalNum) $
    totalNumVec


-- R2S1 both \theta_0 and \theta are represented in the frequency domain
{-# INLINE countR2S1T0T #-}
countR2S1T0T ::
     (Int, Int)
  -> (Int, Int)
  -> [Double]
  -> [Double]
  -> [DList ParticleIndex]
  -> (Int, VU.Vector (Complex Double))
countR2S1T0T (xMin, xMax) (yMin, yMax) t0Freqs tFreqs xs =
  let ys =
        DL.toList .
        DL.concat .
        L.map
          (\((t0f, i), (tf, j)) ->
             DL.map
               (\(x, y, t, _, t0, _) ->
                  ((j, round x, round y, i), exp (0 :+ (t0f * t0 + tf * t)))) .
             DL.concat $
             xs) $
        [(t0f, tf) | t0f <- (L.zip t0Freqs [1 ..]), tf <- (L.zip tFreqs [1 ..])]
      numTrajectories = L.length ys
      arr =
        UA.accum
          (+)
          0
          ((1, xMin, yMin, 1), (L.length tFreqs, xMax, yMax, L.length t0Freqs))
          ys
   in (numTrajectories, toUnboxedVector arr)

-- (Z :. (L.length thetaFreqs) :. xLen :. yLen :. (L.length theta0freqs)) 
solveMonteCarloR2S1T0T ::
     Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> [Double]
  -> [Double]
  -> ParticleIndex
  -> IO (R.Array U DIM4 (Complex Double))
solveMonteCarloR2S1T0T numGen numTrails xLen yLen thetaSigma tao numSteps theta0Freqs thetaFreqs _ = do
  gens <- M.replicateM numGen createSystemRandom
  let !xShift = div xLen 2
      xRange =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      yRange =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 0
  xs <-
    MP.mapM
      (\gen ->
         M.replicateM
           (div numTrails numGen)
           (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
               generatePath
                 gen
                 thetaDist
                 scaleDist
                 10
                 tao
                 numSteps
                 xRange
                 yRange
                 (0, 0, t0, 1, t0, 1)))
      gens
  let (ys, zs) =
        L.unzip $
        parMap rdeepseq (countR2S1T0T xRange yRange theta0Freqs thetaFreqs) xs
      totalNum = fromIntegral $ L.sum ys
      totalNumVec = L.foldl1' (VU.zipWith (+)) zs
  print totalNum
  print $ L.length zs
  return .
    fromUnboxed
      (Z :. (L.length thetaFreqs) :. xLen :. yLen :. (L.length theta0Freqs)) .
    VU.map (\x -> x / totalNum) $
    totalNumVec
