module FokkerPlanck.SpatialSampling
  ( sampleRadial
  , resampleRadial
  ) where

import           Data.Ix
import           Data.List as L

data SpatialSamplingMethod
  = RoundingSS
  | AntialiasingSS
  | GaussianSS Double
               Double
  deriving (Show)

{-# INLINE sampleRadial #-}
sampleRadial :: SpatialSamplingMethod -> (Int, Int) -> Double -> Double -> Bool
sampleRadial RoundingSS rRange x y =
  (round y == 0) && (round x /= 0) && inRange rRange (round x)
sampleRadial AntialiasingSS (rMin, rMax) x y =
  (round x > -1) &&
  inRange (rMin - 1, rMax + 1) (round x) && inRange (-1, 1) (round y)
sampleRadial (GaussianSS radius _) (rMin, rMax) x y =
  (x - (fromIntegral . round $ x)) ^ 2 <= radius ^ 2 - y ^ 2

{-# INLINE distance #-}
distance :: (Double, Double) -> (Double, Double) -> (Double, (Int, Int))
distance (x0, y0) (x1, y1) =
  (sqrt $ (x0 - x1) ^ 2 + (y0 - y1) ^ 2, (round x1, round y1))

{-# INLINE resampleRadial #-}
resampleRadial ::
     SpatialSamplingMethod
  -> (Int, Int)
  -> ((Int, Int, Int, Int, Double, Double), Double)
  -> [((Int, Int, Int, Int, Int), Double)]
resampleRadial RoundingSS _ ((a, b, c, d, x, _), v) = [((a, b, c, d, round x), v)]
resampleRadial AntialiasingSS rRange ((a, b, c, d, x, y), v) =
  let distFunc (x', y') = distance (x, y) (x', y')
      z1@(d1, _) = distFunc (fromIntegral . floor $ x, fromIntegral . floor $ y)
      z2@(d2, _) =
        distFunc (fromIntegral . floor $ x, fromIntegral . ceiling $ y)
      z3@(d3, _) =
        distFunc (fromIntegral . ceiling $ x, fromIntegral . floor $ y)
      z4@(d4, _) =
        distFunc (fromIntegral . ceiling $ x, fromIntegral . ceiling $ y)
      s = d1 + d2 + d3 + d4
  in L.map (\(dist, (x', _)) -> ((a, b, c, d, x'), v * dist / s)) .
     L.filter (\(_, (x', y')) -> y' == 0 && x' /= 0 && inRange rRange x') $
     [z1, z2, z3, z4]
resampleRadial (GaussianSS gRadius sigma) rRange ((a, b, c, d, x, y), v) =
  let deltaX = sqrt $ gRadius ^ 2 - y ^ 2
      xs = [(ceiling $ x - deltaX) .. (floor $ x + deltaX)]
  in L.map
       (\x' ->
          ( (a, b, c, d, x')
          , v * exp (-0.5 * ((x - fromIntegral x') ^ 2 + y ^ 2) / sigma ^ 2))) .
     L.filter (\x' -> x' > 0 && inRange rRange x') $
     xs

{-# INLINE resample #-}
resample ::
     SpatialSamplingMethod
  -> (Int, Int)
  -> (Int, Int)
  -> ((Int, Int, Int, Int, Double, Double), Double)
  -> [((Int, Int, Int, Int, Int, Int), Double)]
resample RoundingSS xRange yRange ((a, b, c, d, x, y), v) =
  let x' = round x
      y' = round y
  in if inRange xRange x' && inRange yRange y' && (x' /= 0 && y' /= 0)
       then [((a, b, c, d, x', y'), v)]
       else []   
resample AntialiasingSS xRange yRange ((a, b, c, d, x, y), v) =
  let distFunc (x', y') = distance (x, y) (x', y')
      z1@(d1, _) = distFunc (fromIntegral . floor $ x, fromIntegral . floor $ y)
      z2@(d2, _) =
        distFunc (fromIntegral . floor $ x, fromIntegral . ceiling $ y)
      z3@(d3, _) =
        distFunc (fromIntegral . ceiling $ x, fromIntegral . floor $ y)
      z4@(d4, _) =
        distFunc (fromIntegral . ceiling $ x, fromIntegral . ceiling $ y)
      s = d1 + d2 + d3 + d4
  in L.map (\(dist, (x', y')) -> ((a, b, c, d, x', y'), v * dist / s)) .
     L.filter
       (\(_, (x', y')) ->
          inRange xRange x' && inRange yRange y' && (x' /= 0 && y' /= 0)) $
     [z1, z2, z3, z4]
resample (GaussianSS gRadius sigma) xRange yRange ((a, b, c, d, x, y), v) =
  let xs =
        [ (x', y')
        | x' <- [(ceiling $ x - gRadius) .. (floor $ x + gRadius)]
        , y' <- [(ceiling $ y - gRadius) .. (floor $ y + gRadius)]
        ]
  in L.map
       (\(x', y') ->
          ( (a, b, c, d, x', y')
          , v *
            exp
              (-0.5 * ((x - fromIntegral x') ^ 2 + (y - fromIntegral y') ^ 2) /
                sigma ^ 2))) .
     L.filter
       (\(x', y') ->
          inRange xRange x' && inRange yRange y' && (x' /= 0 && y' /= 0)) $
     xs
