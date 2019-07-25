module Utils.Coordinates where

{-Coordinate Functions-}
{-# INLINE angleFunctionDeg #-}             

angleFunctionDeg
  :: (Floating a, Ord a)
  => a -> a -> a
angleFunctionDeg 0 0 = 0
angleFunctionDeg i 0
  | i > 0 = 0.0
  | otherwise = 180.0
angleFunctionDeg 0 j
  | j > 0 = 90.0
  | otherwise = 270.0
angleFunctionDeg i j
  | i > 0
  , ratio > 0 = ratio / (pi / 2) * 90.0
  | i > 0
  , ratio < 0 = 360.0 + ratio / (pi / 2) * 90
  | otherwise = 180 + ratio / (pi / 2) * 90
  where
    ratio = atan $ j / i
    

{-# INLINE angleFunctionRad #-}             

angleFunctionRad
  :: (Floating a, Ord a)
  => a -> a -> a
angleFunctionRad 0 0 = 0
angleFunctionRad i 0
  | i > 0 = 0.0
  | otherwise = 1.0 * pi
angleFunctionRad 0 j
  | j > 0 = pi / 2.0
  | otherwise = 3.0 * pi / 2.0
angleFunctionRad i j
  | i > 0
  , ratio > 0 = ratio
  | i > 0
  , ratio < 0 = 2.0 * pi + ratio
  | otherwise = pi + ratio
  where
    ratio = atan $ j / i


{-# INLINE deg2Rad #-}

deg2Rad
  :: (Floating a)
  => a -> a
deg2Rad deg = (deg / 180) * pi

{-# INLINE rad2Deg #-}

rad2Deg
  :: (Floating a)
  => a -> a
rad2Deg rad = (rad / (2 * pi)) * 360.0

coordinateCenter
  :: (Floating a)
  => (a, a) -> (a, a) -> (a, a)
coordinateCenter (i, j) (x, y) = (i - x, j - y)

coordinateCenter'
  :: (Floating a)
  => (a, a) -> (a, a) -> (a, a)
coordinateCenter' (i, j) (x, y) = (i + x, j + y)

cartesian2polar
  :: (Floating a, Ord a)
  => (a, a) -> (a, a)
cartesian2polar (i, j) = (sqrt (i ^ 2 + j ^ 2), angleFunctionDeg i j)

polar2cartesian
  :: (Floating a)
  => (a, a) -> (a, a)
polar2cartesian (rad, deg) =
  (rad * (cos $ deg2Rad deg), rad * (sin $ deg2Rad deg))
  
