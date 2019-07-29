{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.MonteCarlo
  ( module FokkerPlanck.Types
  , module FokkerPlanck.Histogram
  , solveMonteCarloR2S1
  , solveMonteCarloR2S1RP
  , solveMonteCarloR2S1T0
  , solveMonteCarloR2Z1T0
  , solveMonteCarloR2Z2T0S0
  , solveMonteCarloR2Z1T0Radial
  , solveMonteCarloR2Z2T0S0Radial
  , solveMonteCarloR2Z2T0S0Radial'
  , solveMonteCarloR2Z2T0S0ReversalCornerRadial
  ) where

import           Array.UnboxedArray              as UA
import           Control.Arrow
import           Control.Monad                   as M
import           Control.Monad.Parallel          as MP
import           Data.Array                      as Arr
import           Data.Array.Repa                 as R
import           Data.Binary                     (encodeFile)
import           Data.Complex
import           Data.DList                      as DL
import           Data.Ix
import           Data.List                       as L
import           Data.Vector.Unboxed             as VU
import           FokkerPlanck.Histogram
import           FokkerPlanck.Types
import           Statistics.Distribution
import           Statistics.Distribution.Normal
import           Statistics.Distribution.Uniform
import           Statistics.Distribution.Exponential
import           System.Random
import           System.Random.MWC
import           Text.Printf
import           Types
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
    then thetaCheck $ theta + 2 * pi
    else if theta >= 2 * pi
           then thetaCheck $ theta - 2 * pi
           else theta

{-# INLINE scalePlus #-}
scalePlus :: Double -> Double -> Double -> (Double, Bool)
scalePlus maxScale x y -- = (z, False)
  | z < 0 = (0, False)
  -- | z < 0 = (maxScale + z, False)
  | z >= log maxScale = (log maxScale, False)
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
  -> (Int, Int)
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePath randomGen thetaDist' scaleDist' maxScale tao xRange yRange init@(x', y', theta', scale', theta0', scale0') =
  go (x', y', thetaCheck theta', scale', thetaCheck theta0', scale0') DL.empty
  where
    go z@(!x, !y, !theta, !scale, theta0, scale0) xs = do
      deltaTheta <-
        case thetaDist' of
          Nothing -> return 0
          Just thetaDist -> genContVar thetaDist randomGen
      deltaScale <-
        case scaleDist' of
          Nothing -> return 0
          Just scaleDist -> genContVar scaleDist randomGen
      let (!newScale, flag) = scalePlus maxScale scale deltaScale
          eScale = exp scale
          !newTheta =
            if flag
              then (theta `thetaPlus` (deltaTheta + pi))
              else (theta `thetaPlus` deltaTheta)
          !newX = x + eScale * cos theta
          !newY = y + eScale * sin theta
          newIndex = (newX, newY, newTheta, newScale, theta0, scale0)
          ys =
            if (inRange xRange (round newX)) && (inRange yRange (round newY))
              then DL.cons newIndex xs
              else xs
      t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
      if t < (1 - exp ((-1) / tao))
        then return ys
        else go newIndex ys

-- R2S1

{-# INLINE countR2S1 #-}
countR2S1 ::
     (Int, Int)
  -> (Int, Int)
  -> Int
  -> [DList ParticleIndex]
  -> Histogram Int
countR2S1 (xMin, xMax) (yMin, yMax) numOrientations xs =
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        L.map (\(t, x, y) -> ((floor $ t / deltaTheta, round x, round y), 1)) .
        DL.toList .
        DL.concat . L.map (DL.map (\(x, y, t, s, t0, s0) -> (t, x, y))) $
        xs
      numTrajectories = L.length ys
      arr =
        UA.accum (+) 0 ((0, xMin, yMin), (numOrientations - 1, xMax, yMax)) ys
   in Histogram
        [(yMax - yMin + 1), (xMax - xMin + 1), numOrientations]
        numTrajectories
        (toUnboxedVector arr)

solveMonteCarloR2S1 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> FilePath
  -> ParticleIndex
  -> IO R2S1Array
solveMonteCarloR2S1 numGen numTrails xLen yLen numOrientations thetaSigma tao r histFilePath (_, _, _, _, t0, s0) = do
  gens <- M.replicateM numGen createSystemRandom
  let !xShift = div xLen 2
      func (a, b) =
        ( if abs a > r
            then -r
            else a
        , if b > r
            then r
            else b)
      xRange' =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      !yShift = div yLen 2
      yRange' =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      xRange = func xRange'
      yRange = func yRange'
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 0
  xs <-
    MP.mapM
      (\gen ->
         M.replicateM
           (div numTrails numGen)
           (do deltaTheta <-
                 case thetaDist of
                   Nothing -> return 0
                   Just thetaDist' -> genContVar thetaDist' gen
               generatePath
                 gen
                 thetaDist
                 scaleDist
                 10
                 tao
                 xRange
                 yRange
                 (0, 0, t0 + deltaTheta, s0, t0 + deltaTheta, s0)))
      gens
  let ys = parMap rdeepseq (countR2S1 xRange' yRange' numOrientations) xs
      histogram = L.foldl1' addHistogram ys
  unless (L.null histFilePath) (encodeFile histFilePath histogram)
  return . getNormalizedHistogramArr . mapHistogram (\x -> fromIntegral x :+ 0) $
    histogram

-- R2S1RP

{-# INLINE countR2S1RP #-}
countR2S1RP ::
     (Int, Int)
  -> (Int, Int)
  -> Int
  -> Int
  -> Double
  -> [DList ParticleIndex]
  -> Histogram Int
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
   in Histogram
        [(yMax - yMin + 1), (xMax - xMin + 1), numScales, numOrientations]
        numTrajectories
        (toUnboxedVector arr)

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
  -> IO R2S1RPArray
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
              xRange
              yRange
              (0, 0, t0, s0, t0, s0)))
      gens
  let ys =
        parMap
          rdeepseq
          (countR2S1RP xRange yRange numOrientations numScales maxScale)
          xs
      histogram = L.foldl1' addHistogram ys
  return . getNormalizedHistogramArr . mapHistogram fromIntegral $ histogram

-- R2S1 \theta_0 is represented in the frequency domain
{-# INLINE countR2S1T0 #-}
countR2S1T0 ::
     (Int, Int)
  -> (Int, Int)
  -> Int
  -> [Double]
  -> [DList ParticleIndex]
  -> Histogram (Complex Double)
countR2S1T0 (xMin, xMax) (yMin, yMax) numOrientations freqs xs =
  let !deltaTheta = 2 * pi / (fromIntegral numOrientations)
      ys =
        DL.toList . DL.concat $
        L.zipWith
          (\freq i ->
             DL.map
               (\(x, y, t, _, t0, _) ->
                  ( (i, floor (t / deltaTheta), round x, round y)
                  , exp (0 :+ (freq * t0)))) .
             DL.concat $
             xs)
          freqs
          [1 ..]
      numTrajectories = div (L.length ys) (L.length freqs)
      arr =
        UA.accum
          (+)
          0
          ( (1, 0, xMin, yMin)
          , (L.length freqs, numOrientations - 1, xMax, yMax))
          ys
   in Histogram
        [ (yMax - yMin + 1)
        , (xMax - xMin + 1)
        , numOrientations
        , (L.length freqs)
        ]
        numTrajectories
        (toUnboxedVector arr)

-- (Z :. (L.length freqs) :. numOrientations  :. xLen :. yLen )
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
  -> IO R2S1T0Array
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
                 xRange
                 yRange
                 (0, 0, t0, 1, t0, 1)))
      gens
  let ys =
        parMap
          rdeepseq
          (countR2S1T0 xRange yRange numOrientations theta0freqs)
          xs
      histogram = L.foldl1' addHistogram ys
  return . getNormalizedHistogramArr $ histogram

-- R2S1 both \theta_0 and \theta are represented in the frequency domain
{-# INLINE countR2Z1T0 #-}
countR2Z1T0 ::
     (Int, Int)
  -> (Int, Int)
  -> [Double]
  -> [Double]
  -> [DList ParticleIndex]
  -> Histogram (Complex Double)
countR2Z1T0 (xMin, xMax) (yMin, yMax) t0Freqs tFreqs xs =
  let ys =
        DL.toList .
        DL.concat .
        L.map
          (\((t0f, i), (tf, j)) ->
             DL.map
               (\(x, y, t, _, t0, _) ->
                  ((j, i, round x, round y), exp (0 :+ (t0f * t0 + tf * t)))) .
             DL.concat $
             xs) $
        [(t0f, tf) | t0f <- (L.zip t0Freqs [1 ..]), tf <- (L.zip tFreqs [1 ..])]
      numTrajectories =
        div (L.length ys) ((L.length t0Freqs) * (L.length tFreqs))
      arr =
        UA.accum
          (+)
          0
          ((1, 1, xMin, yMin), (L.length tFreqs, L.length t0Freqs, xMax, yMax))
          ys
   in Histogram
        [ (yMax - yMin + 1)
        , (xMax - xMin + 1)
        , (L.length t0Freqs)
        , (L.length tFreqs)
        ]
        numTrajectories
        (toUnboxedVector arr)

-- (Z :. (L.length thetafreqs) :. (L.length theta0Freqs) :. xLen :. yLen )
solveMonteCarloR2Z1T0 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Int
  -> [Double]
  -> [Double]
  -> FilePath
  -> Histogram (Complex Double)
  -> IO R2Z1T0Array
solveMonteCarloR2Z1T0 numGen numTrails maxTrails xLen yLen thetaSigma tao initialScale numSteps theta0Freqs thetaFreqs filePath hist = do
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
      numMonteCarlo = div numTrails maxTrails
      numLeft = mod numTrails maxTrails
  histograms <-
    M.replicateM
      numMonteCarlo
      (do xs <-
            MP.mapM
              (\gen ->
                 M.replicateM
                   (div maxTrails numGen)
                   (do t0 <-
                         genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                       generatePath
                         gen
                         thetaDist
                         scaleDist
                         initialScale
                         tao
                         xRange
                         yRange
                         (0, 0, t0, initialScale, t0, initialScale)))
              gens
          let ys =
                parMap
                  rdeepseq
                  (countR2Z1T0 xRange yRange theta0Freqs thetaFreqs)
                  xs
              !histogram = L.foldl1' addHistogram ys
          return $! histogram)
  if numLeft > 0
    then do
      xs <-
        MP.mapM
          (\gen ->
             M.replicateM
               (div numLeft numGen)
               (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                   generatePath
                     gen
                     thetaDist
                     scaleDist
                     initialScale
                     tao
                     xRange
                     yRange
                     (0, 0, t0, initialScale, t0, initialScale)))
          gens
      let ys =
            parMap
              rdeepseq
              (countR2Z1T0 xRange yRange theta0Freqs thetaFreqs)
              xs
          histogram' = L.foldl1' addHistogram ys
          histogram = L.foldl' addHistogram histogram' histograms
          !newHist = addHistogram histogram hist
      unless (L.null filePath) (encodeFile filePath newHist)
      return . getNormalizedHistogramArr $ newHist
    else do
      let !histogram = L.foldl' addHistogram hist histograms
      unless (L.null filePath) (encodeFile filePath histogram)
      return . getNormalizedHistogramArr $ histogram


-- R2S1RP \theta_0, \theta, s0 and s are represented in the frequency domain
{-# INLINE countR2Z2T0S0 #-}
countR2Z2T0S0 ::
     (Int, Int)
  -> (Int, Int)
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> [DList ParticleIndex]
  -> Histogram (Complex Double)
countR2Z2T0S0 (xMin, xMax) (yMin, yMax) t0Freqs tFreqs s0Freqs sFreqs maxR xs =
  let ys =
        DL.toList .
        DL.concat .
        L.map
          (\((t0f, i), (tf, j), (s0f, k), (sf, l)) ->
             DL.map
               (\(x, y, t, s, t0, s0) ->
                  let -- s' =
                      --   if s == 0
                      --     then 0
                      --     else log s
                      -- s0' =
                      --   if s0 == 0
                      --     then 0
                      --     else log s0
                      !v =
                        exp
                          (0 :+
                           (-t0f * t0 + tf * t +
                            (sf * s + s0f * s0) * 2 * pi / log maxR))
                      !x' = round x
                      !y' = round y
                   in ((j, l, i, k, x', y'), v)) .
             DL.concat $
             xs) $
        [ (t0f, tf, s0f, sf)
        | t0f <- (L.zip t0Freqs [1 ..])
        , tf <- (L.zip tFreqs [1 ..])
        , s0f <- (L.zip s0Freqs [1 ..])
        , sf <- (L.zip sFreqs [1 ..])
        ]
      numTrajectories =
        div
          (L.length ys)
          ((L.length t0Freqs) * (L.length tFreqs) * (L.length s0Freqs) *
           (L.length sFreqs))
      arr =
        UA.accum
          (+)
          0
          ( (1, 1, 1, 1, xMin, yMin)
          , ( L.length tFreqs
            , L.length sFreqs
            , L.length t0Freqs
            , L.length s0Freqs
            , xMax
            , yMax))
          ys
   in Histogram
        [ (yMax - yMin + 1)
        , (xMax - xMin + 1)
        , (L.length s0Freqs)
        , (L.length t0Freqs)
        , (L.length sFreqs)
        , (L.length tFreqs)
        ]
        numTrajectories
        (toUnboxedVector arr)


-- (Z :. (L.length thetafreqs) :. (L.length scaleFreqs) :. (L.length theta0Freqs) :. (L.length scale0Freqs) :. xLen :. yLen )
solveMonteCarloR2Z2T0S0 ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> FilePath
  -> Histogram (Complex Double)
  -> IO R2Z2T0S0Array
solveMonteCarloR2Z2T0S0 numGen numTrails maxTrails xLen yLen thetaSigma scaleSigma maxScale tao numSteps theta0Freqs thetaFreqs scale0Freqs scaleFreqs filePath hist = do
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
      numMonteCarlo = div numTrails maxTrails
      numLeft = mod numTrails maxTrails
  printf
    "%d trails, maximum %d each batch, %d batch in total.\n"
    numTrails
    maxTrails
    (numMonteCarlo +
     if numLeft > 0
       then 1
       else 0)
  gensList <-
    M.replicateM numMonteCarlo (M.replicateM numGen createSystemRandom)
  histogram <-
    M.foldM
      (\h gens -> do
         xs <-
           MP.mapM
             (\gen ->
                M.replicateM
                  (div maxTrails numGen)
                  (do t0 <-
                        genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                      s0 <-
                        genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                      generatePath
                        gen
                        thetaDist
                        scaleDist
                        maxScale
                        tao
                        xRange
                        yRange
                        (0, 0, t0, s0, t0, s0)))
             gens
         let ys =
               parMap
                 rdeepseq
                 (countR2Z2T0S0
                    xRange
                    yRange
                    theta0Freqs
                    thetaFreqs
                    scale0Freqs
                    scaleFreqs
                    maxScale)
                 xs
         return $! L.foldl' addHistogram h ys)
      hist
      gensList
  if numLeft > 0
    then do
      gens <- M.replicateM numGen createSystemRandom
      xs <-
        MP.mapM
          (\gen ->
             M.replicateM
               (div numLeft numGen)
               (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                   s0 <-
                     genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                   generatePath
                     gen
                     thetaDist
                     scaleDist
                     maxScale
                     tao
                     xRange
                     yRange
                     (0, 0, t0, s0, t0, s0)))
          gens
      let ys =
            parMap
              rdeepseq
              (countR2Z2T0S0
                 xRange
                 yRange
                 theta0Freqs
                 thetaFreqs
                 scale0Freqs
                 scaleFreqs
                 maxScale)
              xs
          histogram' = L.foldl' addHistogram histogram ys
          newHist =
            if numMonteCarlo > 0
              then histogram'
              else addHistogram hist histogram'
      unless (L.null filePath) (encodeFile filePath newHist)
      return . getNormalizedHistogramArr $ newHist
    else do
      unless (L.null filePath) (encodeFile filePath histogram)
      return . getNormalizedHistogramArr $ histogram




-- Only simulate the radial part
{-# INLINE generatePathRaidal #-}
generatePathRaidal ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePathRaidal randomGen thetaDist' scaleDist' maxScale tao rRange init@(x', y', theta', scale', theta0', scale0') =
  go (x', y', thetaCheck theta', scale', thetaCheck theta0', scale0') DL.empty
  where
    go z@(!x, !y, !theta, !scale, theta0, scale0) xs = do
      deltaTheta <-
        case thetaDist' of
          Nothing -> return 0
          Just thetaDist -> genContVar thetaDist randomGen
      deltaScale <-
        case scaleDist' of
          Nothing -> return 0
          Just scaleDist -> genContVar scaleDist randomGen
      let (newScale, _) = scalePlus maxScale scale deltaScale
          newTheta = theta `thetaPlus` deltaTheta
          newX = x + (exp scale) * cos theta
          newY = y + (exp scale) * sin theta
          -- r = sqrt $ newX' ^ 2 + newY' ^ 2
          -- t = angleFunctionDeg newX' newY'
          -- (newX, newY) =
          --   if r > maxScale
          --     then ( (2 * maxScale - r) * sin (t + pi)
          --          , (2 * maxScale - r) * cos (t + pi))
          --     else (newX', newY')
          newIndex = (newX, newY, newTheta, newScale, theta0, scale0)
          recordIndex = (newX, newY, theta, scale, theta0, scale0)
          ys =
            if (inRange rRange (round newX)) && round newY == 0
              then DL.cons newIndex xs
              else xs
          -- ys =
          --   if (inRange rRange (round x)) && round y == 0
          --     then DL.cons z xs
          --     else xs
      t <- genContVar (uniformDistr 0 1) randomGen :: IO Double
      if t > exp ((-1) / tao)
        then return ys
        else go newIndex ys


-- Raidal version: R2S1 both \theta_0 and \theta are represented in the frequency domain
{-# INLINE countR2Z1T0Radial #-}
countR2Z1T0Radial ::
     (Int, Int)
  -> [Double]
  -> [Double]
  -> [DList ParticleIndex]
  -> Histogram (Complex Double)
countR2Z1T0Radial (rMin, rMax) t0Freqs tFreqs xs =
  let ys =
        DL.toList .
        DL.concat .
        L.map
          (\((t0f, i), (tf, j)) ->
             DL.map
               (\(x, _, t, _, t0, _) ->
                  ((j, i, round x), exp (0 :+ (-t0f * t0 + tf * t)))) .
             DL.concat $
             xs) $
        [(t0f, tf) | t0f <- (L.zip t0Freqs [1 ..]), tf <- (L.zip tFreqs [1 ..])]
      numTrajectories =
        div (L.length ys) ((L.length t0Freqs) * (L.length tFreqs))
      arr =
        UA.accum
          (+)
          0
          ((1, 1, rMin), (L.length tFreqs, L.length t0Freqs, rMax))
          ys
   in Histogram
        [(rMax - rMin + 1), (L.length t0Freqs), (L.length tFreqs)]
        numTrajectories
        (toUnboxedVector arr)

-- Radial vestion: (Z :. (L.length thetafreqs) :. (L.length theta0Freqs) :. xLen)
solveMonteCarloR2Z1T0Radial ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (R.Array D DIM3 Double)
solveMonteCarloR2Z1T0Radial numGen numTrails maxTrails xLen yLen thetaSigma tao initialScale theta0Freqs thetaFreqs filePath hist = do
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
      rRange = (0, (round . sqrt . fromIntegral $ xShift ^ 2 + yShift ^ 2) - 1)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 0
      numMonteCarlo = div numTrails maxTrails
      numLeft = mod numTrails maxTrails
  histograms <-
    M.replicateM
      numMonteCarlo
      (do xs <-
            MP.mapM
              (\gen ->
                 M.replicateM
                   (div maxTrails numGen)
                   (do t0 <-
                         genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                       generatePathRaidal
                         gen
                         thetaDist
                         scaleDist
                         initialScale
                         tao
                         rRange
                         (0, 0, t0, initialScale, t0, initialScale)))
              gens
          let ys =
                parMap
                  rdeepseq
                  (countR2Z1T0Radial rRange theta0Freqs thetaFreqs)
                  xs
              !histogram = L.foldl1' addHistogram ys
          return $! histogram)
  if numLeft > 0
    then do
      xs <-
        MP.mapM
          (\gen ->
             M.replicateM
               (div numLeft numGen)
               (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                   generatePathRaidal
                     gen
                     thetaDist
                     scaleDist
                     initialScale
                     tao
                     rRange
                     (0, 0, t0, initialScale, t0, initialScale)))
          gens
      let ys =
            parMap rdeepseq (countR2Z1T0Radial rRange theta0Freqs thetaFreqs) xs
          histogram' = L.foldl1' addHistogram ys
          histogram = L.foldl' addHistogram histogram' histograms
          !newHist = addHistogram histogram hist
      unless (L.null filePath) (encodeFile filePath newHist)
      return . R.map magnitude . getNormalizedHistogramArr $ newHist
    else do
      let !histogram = L.foldl' addHistogram hist histograms
      unless (L.null filePath) (encodeFile filePath histogram)
      return . R.map magnitude . getNormalizedHistogramArr $ histogram


-- Radial version: R2S1RP \theta_0, \theta, s0 and s are represented in the frequency domain
{-# INLINE countR2Z2T0S0Radial #-}
countR2Z2T0S0Radial ::
     (Int, Int)
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> [DList ParticleIndex]
  -> Histogram (Complex Double)
countR2Z2T0S0Radial (rMin, rMax) t0Freqs tFreqs s0Freqs sFreqs maxScale xs =
  let ys =
        DL.toList .
        DL.concat .
        L.map
          (\((t0f, i), (tf, j), (s0f, k), (sf, l)) ->
             DL.map
               (\(x, _, t, s, t0, s0) ->
                  let -- s' =
                      --   if s == 0
                      --     then 0
                      --     else log s
                      -- s0' =
                      --   if s0 == 0
                      --     then 0
                      --     else log s0
                      !v =
                        exp $
                        0 :+
                        (-t0f * t0 + tf * t +
                         (sf * s + s0f * s0) * 2 * pi / (log maxScale))
                      !x' = round x
                   in ((j, l, i, k, x'), v)) .
             DL.concat $
             xs) $
        [ (t0f, tf, s0f, sf)
        | t0f <- (L.zip t0Freqs [1 ..])
        , tf <- (L.zip tFreqs [1 ..])
        , s0f <- (L.zip s0Freqs [1 ..])
        , sf <- (L.zip sFreqs [1 ..])
        ]
      numTrajectories =
        div
          (L.length ys)
          ((L.length t0Freqs) * (L.length tFreqs) * (L.length s0Freqs) *
           (L.length sFreqs))
      arr =
        UA.accum
          (+)
          0
          ( (1, 1, 1, 1, rMin)
          , ( L.length tFreqs
            , L.length sFreqs
            , L.length t0Freqs
            , L.length s0Freqs
            , rMax))
          ys
   in Histogram
        [ (rMax - rMin + 1)
        , (L.length s0Freqs)
        , (L.length t0Freqs)
        , (L.length sFreqs)
        , (L.length tFreqs)
        ]
        numTrajectories
        (toUnboxedVector arr)


-- (Z :. (L.length thetafreqs) :. (L.length scaleFreqs) :. (L.length theta0Freqs) :. (L.length scale0Freqs) :. rLen )
solveMonteCarloR2Z2T0S0Radial ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (R.Array D DIM5 Double)
solveMonteCarloR2Z2T0S0Radial numGen numTrails maxTrails xLen yLen thetaSigma scaleSigma maxScale tao theta0Freqs thetaFreqs scale0Freqs scaleFreqs filePath hist = do
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
      rRange = (0, (round . sqrt . fromIntegral $ xShift ^ 2 + yShift ^ 2) - 1)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 scaleSigma
      numMonteCarlo = div numTrails maxTrails
      numLeft = mod numTrails maxTrails
  printf
    "%d trails, maximum %d each batch, %d batch in total.\n"
    numTrails
    maxTrails
    (numMonteCarlo +
     if numLeft > 0
       then 1
       else 0)
  gensList <-
    M.replicateM numMonteCarlo (M.replicateM numGen createSystemRandom)
  histogram <-
    M.foldM
      (\h gens -> do
         xs <-
           MP.mapM
             (\gen ->
                M.replicateM
                  (div maxTrails numGen)
                  (do t0 <-
                        genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                      s0 <-
                        genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                      generatePathRaidal
                        gen
                        thetaDist
                        scaleDist
                        maxScale
                        tao
                        rRange
                        (0, 0, t0, s0, t0, s0)))
             gens 
         let ys =
               parMap
                 rdeepseq
                 (countR2Z2T0S0Radial
                    rRange
                    theta0Freqs
                    thetaFreqs
                    scale0Freqs
                    scaleFreqs
                    maxScale)
                 xs
         return $! L.foldl' addHistogram h ys)
      hist
      gensList
  if numLeft > 0
    then do
      gens <- M.replicateM numGen createSystemRandom
      xs <-
        MP.mapM
          (\gen ->
             M.replicateM
               (div numLeft numGen)
               (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                   s0 <-
                     genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                   generatePathRaidal
                     gen
                     thetaDist
                     scaleDist
                     maxScale
                     tao
                     rRange
                     (0, 0, t0, s0, t0, s0)))
          gens
      let ys =
            parMap
              rdeepseq
              (countR2Z2T0S0Radial
                 rRange
                 theta0Freqs
                 thetaFreqs
                 scale0Freqs
                 scaleFreqs
                 maxScale
                 )
              xs
          histogram' = L.foldl' addHistogram histogram ys
          newHist =
            if numMonteCarlo > 0
              then histogram'
              else addHistogram hist histogram'
      unless (L.null filePath) (encodeFile filePath newHist)
      return . R.map magnitude . getNormalizedHistogramArr $ newHist
    else do
      unless (L.null filePath) (encodeFile filePath histogram)
      return . R.map magnitude . getNormalizedHistogramArr $ histogram
      

solveMonteCarloR2Z2T0S0Radial' ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (R.Array U DIM5 (Complex Double))
solveMonteCarloR2Z2T0S0Radial' numGen numTrails maxTrails xLen yLen thetaSigma scaleSigma maxScale tao theta0Freqs thetaFreqs scale0Freqs scaleFreqs filePath hist = do
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
      rRange = (0, (round . sqrt . fromIntegral $ xShift ^ 2 + yShift ^ 2) - 1)
      -- rRange = (0, round maxScale)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 scaleSigma
      numMonteCarlo = div numTrails maxTrails
      numLeft = mod numTrails maxTrails
  printf
    "%d trails, maximum %d each batch, %d batch in total.\n"
    numTrails
    maxTrails
    (numMonteCarlo +
     if numLeft > 0
       then 1
       else 0)
  gensList <-
    M.replicateM numMonteCarlo (M.replicateM numGen createSystemRandom)
  histogram <-
    M.foldM
      (\h gens -> do
         xs <-
           MP.mapM
             (\gen ->
                M.replicateM
                  (div maxTrails numGen)
                  (do t0 <-
                        genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                      s0 <-
                        genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                      generatePathRaidal
                        gen
                        thetaDist
                        scaleDist
                        maxScale
                        tao
                        rRange
                        (0, 0, t0, s0, t0, s0)))
             gens
         let ys =
               parMap
                 rdeepseq
                 (countR2Z2T0S0Radial
                    rRange
                    theta0Freqs
                    thetaFreqs
                    scale0Freqs
                    scaleFreqs
                    maxScale)
                 xs
         return $! L.foldl' addHistogram h ys)
      hist
      gensList
  if numLeft > 0
    then do
      gens <- M.replicateM numGen createSystemRandom
      xs <-
        MP.mapM
          (\gen ->
             M.replicateM
               (div numLeft numGen)
               (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                   s0 <-
                     genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                   generatePathRaidal
                     gen
                     thetaDist
                     scaleDist
                     maxScale
                     tao
                     rRange
                     (0, 0, t0, s0, t0, s0)))
          gens
      let ys =
            parMap
              rdeepseq
              (countR2Z2T0S0Radial
                 rRange
                 theta0Freqs
                 thetaFreqs
                 scale0Freqs
                 scaleFreqs
                 maxScale)
              xs
          histogram' = L.foldl' addHistogram histogram ys
          newHist =
            if numMonteCarlo > 0
              then histogram'
              else addHistogram hist histogram'
      -- unless (L.null filePath) (encodeFile filePath newHist)
      return . getNormalizedHistogramArrSink' $ newHist
    else do
      -- unless (L.null filePath) (encodeFile filePath histogram)
      return . getNormalizedHistogramArrSink' $ histogram




{-# INLINE generatePathReversalCornerRaidal #-}
generatePathReversalCornerRaidal ::
     (Distribution d, ContGen d)
  => GenIO
  -> Maybe d
  -> Maybe d
  -> Double
  -> Double
  -> Double
  -> Double
  -> (Int, Int)
  -> ParticleIndex
  -> IO (DList ParticleIndex)
generatePathReversalCornerRaidal randomGen thetaDist' scaleDist' maxScale taoDecay taoReversal taoCorner rRange init@(x', y', theta', scale', theta0', scale0') =
  go (x', y', thetaCheck theta', scale', thetaCheck theta0', scale0') DL.empty
  where
    go z@(!x, !y, !theta, !scale, theta0, scale0) xs = do
      deltaTheta <-
        case thetaDist' of
          Nothing -> return 0
          Just thetaDist -> genContVar thetaDist randomGen
      deltaScale <-
        case scaleDist' of
          Nothing -> return 0
          Just scaleDist -> genContVar scaleDist randomGen
      let newScale = scale + deltaScale
          newTheta' = theta `thetaPlus` deltaTheta
          newX = x + (exp scale) * cos theta
          newY = y + (exp scale) * sin theta
      newTheta <-
        do tReversal <- genContVar (uniformDistr 0 1) randomGen :: IO Double
           if tReversal < 1 - exp ((-1) / taoReversal)
             then return $ newTheta' `thetaPlus` pi
             else return newTheta' >>= \t' -> do
                    tCorner <-
                      genContVar (uniformDistr 0 1) randomGen :: IO Double
                    if tCorner < 1 - exp ((-1) / taoCorner)
                      then thetaCheck <$>
                           genContVar (uniformDistr 0 (2 * pi)) randomGen :: IO Double
                      else return t'
      let newIndex = (newX, newY, newTheta, newScale, theta0, scale0)
          ys =
            if (inRange rRange (round newX)) && round newY == 0
              then DL.cons newIndex xs
              else xs
      tDecay <- genContVar (uniformDistr 0 1) randomGen :: IO Double
      if tDecay < 1 - exp ((-1) / taoDecay)
        then return ys
        else go newIndex ys 

solveMonteCarloR2Z2T0S0ReversalCornerRadial ::
     Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> FilePath
  -> Histogram (Complex Double)
  -> IO (R.Array D DIM5 Double)
solveMonteCarloR2Z2T0S0ReversalCornerRadial numGen numTrails maxTrails xLen yLen thetaSigma scaleSigma maxScale taoDecay taoReversal taoCorner theta0Freqs thetaFreqs scale0Freqs scaleFreqs filePath hist = do
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
      rRange = (0, (round . sqrt . fromIntegral $ xShift ^ 2 + yShift ^ 2) - 1)
      thetaDist = normalDistrE 0 thetaSigma
      scaleDist = normalDistrE 0 scaleSigma
      numMonteCarlo = div numTrails maxTrails
      numLeft = mod numTrails maxTrails
  printf
    "%d trails, maximum %d each batch, %d batch in total.\n"
    numTrails
    maxTrails
    (numMonteCarlo +
     if numLeft > 0
       then 1
       else 0)
  gensList <-
    M.replicateM numMonteCarlo (M.replicateM numGen createSystemRandom)
  histogram <-
    M.foldM
      (\h gens -> do
         xs <-
           MP.mapM
             (\gen ->
                M.replicateM
                  (div maxTrails numGen)
                  (do t0 <-
                        genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                      s0 <-
                        genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                      generatePathReversalCornerRaidal
                        gen
                        thetaDist
                        scaleDist
                        maxScale
                        taoDecay
                        taoReversal
                        taoCorner
                        rRange
                        (0, 0, t0, s0, t0, s0)))
             gens
         let ys =
               parMap
                 rdeepseq
                 (countR2Z2T0S0Radial
                    rRange
                    theta0Freqs
                    thetaFreqs
                    scale0Freqs
                    scaleFreqs
                    maxScale)
                 xs
         return $! L.foldl' addHistogram h ys)
      hist
      gensList
  if numLeft > 0
    then do
      gens <- M.replicateM numGen createSystemRandom
      xs <-
        MP.mapM
          (\gen ->
             M.replicateM
               (div numLeft numGen)
               (do t0 <- genContVar (uniformDistr 0 (2 * pi)) gen :: IO Double
                   s0 <-
                     genContVar (uniformDistr 0 (log maxScale)) gen :: IO Double
                   generatePathReversalCornerRaidal
                     gen
                     thetaDist
                     scaleDist
                     maxScale
                     taoDecay
                     taoReversal
                     taoCorner
                     rRange
                     (0, 0, t0, s0, t0, s0)))
          gens
      let ys =
            parMap
              rdeepseq
              (countR2Z2T0S0Radial
                 rRange
                 theta0Freqs
                 thetaFreqs
                 scale0Freqs
                 scaleFreqs
                 maxScale)
              xs
          histogram' = L.foldl' addHistogram histogram ys
          newHist =
            if numMonteCarlo > 0
              then histogram'
              else addHistogram hist histogram'
      unless (L.null filePath) (encodeFile filePath newHist)
      return . R.map magnitude . getNormalizedHistogramArr $ newHist
    else do
      unless (L.null filePath) (encodeFile filePath histogram)
      return . R.map magnitude . getNormalizedHistogramArr $ histogram
