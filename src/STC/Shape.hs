module STC.Shape where

import           Data.Array.Repa    as R
import           Data.Array.Unboxed as AU
import           Data.Ix
import           Data.List          as L
-- import           Debug.Trace
-- import           Text.Printf
import           Utils.Coordinates

data Points a =
  Points (Double, Double)
         Int
         a
  deriving (Show, Read)

data ContrastArea a =
  ContrastArea Bool
               a
  deriving (Show)

data Shape2D
  = Line { lineOrientationDeg :: Double
         , lineLength         :: Int
         , lineWidth          :: Int }
  | Corner { cornerThetaDeg0 :: Double
           , cornerThetaDeg1 :: Double
           , cornerLength    :: Int }
  | IncompleteCircle { iCircleTheta0 :: Double
                     , iCircleTheta1 :: Double
                     , iCircleRadiau :: Int }
  | ETriangle { eTriangleThetaDeg :: Double
              , eTriangleLength   :: Int }
  | PacMan { pacManThetaDeg0 :: Double
           , pacManThetaDeg1 :: Double
           , pacManRadius    :: Int }
  | KanizsaTriangle1 { kanizsaTriangle1ThetaDeg        :: Double
                     , kanizsaTriangle1TriangleLength0 :: Int
                     , kanizsaTriangle1TriangleLength1 :: Int
                     , kanizsaTriangle1Radius          :: Int }
  | TJunction { tJunctionThetaDeg :: Double
              , tJunctionLength   :: Int }
  | Cross { crossThetaDeg :: Double
          , crossLength   :: Int }
  | Rectangle { rectangleThetaDeg :: Double
              , rectangleWidth    :: Int
              , rectangleHeight   :: Int }
  | Ehrenstein { ehrensteinNumPoint    :: Int
               , ehrensteinInnerRadius :: Int
               , ehrensteinLength      :: Int }
  | Circle { circleNum    :: Int
           , circleRadius :: Int }
  deriving (Show, Read)

isDuplicate :: Int -> (Double,Double) -> [(Double,Double)] -> Bool
isDuplicate len (x1, y1) =
  L.any
    (\(x2, y2) ->
       let r = sqrt $ (x1 - x2) ^ 2 + (y1 - y2) ^ 2
        in r < (fromIntegral len - 1))

{-# INLINE removeDuplicate #-}
removeDuplicate :: Int -> [(Double, Double)] -> [(Double, Double)]
removeDuplicate _ [] = []
removeDuplicate len (x:xs) =
  if isDuplicate len x xs
    then removeDuplicate len xs
    else x : removeDuplicate len xs

{-# INLINE getShape2DIndexList #-}
getShape2DIndexList :: [(Double, Double)] -> [(Int, Int)]
getShape2DIndexList  =
  L.map (\(x, y) -> (round x, round y))

-- {-# INLINE getShape2DIndexListDouble #-}
-- getShape2DIndexListDouble :: [(Double, Double)] -> [(Double, Double)]
-- getShape2DIndexListDouble = removeDuplicate

{-# INLINE getShape2DRepaArray #-}
getShape2DRepaArray ::
     Int -> Int -> [(Double, Double)] -> R.Array U DIM3 Double
getShape2DRepaArray rows cols xs =
  fromListUnboxed (Z :. (1 :: Int) :. cols :. rows) . AU.elems $
  (accumArray (+) 0 ((0, 0), (cols - 1, rows - 1)) .
   L.map (\x -> (x, 1)) .
   L.filter (\(x, y) -> (inRange (0, cols - 1) x) && (inRange (0, rows - 1) y)) .
   getShape2DIndexList $
   xs :: UArray (Int, Int) Double)

{-# INLINE makeShape2DList #-}
makeShape2DList :: [Points Shape2D] -> [(Double,Double)]
makeShape2DList = L.concatMap makeShape2D

{-# INLINE makeShape2D #-}
makeShape2D :: Points Shape2D -> [(Double,Double)]
makeShape2D (Points (centerC, centerR) minDist (Line ori len _)) =
  removeDuplicate minDist . L.map (\(x, y) -> (x + centerC, y + centerR)) $
  [ (r * cos (deg2Rad ori), r * sin (deg2Rad ori))
  | r <- (L.map fromIntegral [0,minDist .. len])
  ]
makeShape2D (Points (centerC, centerR) minDist (Corner thetaDeg0 thetaDeg1 len)) =
  let thetaRad0 = deg2Rad thetaDeg0
      thetaRad1 = deg2Rad thetaDeg1
      pointLens = L.map fromIntegral [0,minDist .. len]
   in removeDuplicate minDist . L.map (\(x, y) -> (x + centerC, y + centerR)) $
      [ ( r * cos (thetaRad0 - thetaRad1 / 2)
        , r * sin (thetaRad0 - thetaRad1 / 2))
      | r <- pointLens
      ] L.++
      [ ( r * cos (thetaRad0 + thetaRad1 / 2)
        , r * sin (thetaRad0 + thetaRad1 / 2))
      | r <- pointLens
      ]
makeShape2D (Points (centerC, centerR) minDist (ETriangle thetaDeg len)) =
  let x = fromIntegral len / 2 / cos (pi / 6)
      theta = deg2Rad (240 + thetaDeg)
      thetaRad = deg2Rad thetaDeg
   in removeDuplicate minDist . makeShape2DList $
      [ (Points
           (x * cos thetaRad, x * sin thetaRad)
           minDist
           (Corner (180 + thetaDeg) 60 len))
      , (Points
           (x * cos theta, x * sin theta)
           minDist
           (Line (90 + thetaDeg) len 0))
      ]
makeShape2D (Points (centerC, centerR) minDist (IncompleteCircle thetaDeg0 thetaDeg1 r)) =
  let thetaRad0 = deg2Rad thetaDeg0
      thetaRad1 = deg2Rad thetaDeg1
      n =
        floor ((2 * pi - thetaRad1) * fromIntegral r / fromIntegral minDist) + 1 :: Int
      deltaTheta = (2 * pi - thetaRad1) / fromIntegral n
      xs =
        [ (fromIntegral r * cos theta, fromIntegral r * sin theta)
        | theta <-
            [(thetaRad0 + thetaRad1 / 2),(thetaRad0 + thetaRad1 / 2 + deltaTheta) .. (thetaRad0 -
                                                                                      thetaRad1 /
                                                                                      2 +
                                                                                      2 *
                                                                                      pi)]
        ]
   in removeDuplicate minDist . L.map (\(x, y) -> (x + centerC, y + centerR)) $
      xs
makeShape2D (Points (centerC, centerR) minDist (PacMan thetaDeg0 thetaDeg1 r)) =
  removeDuplicate minDist . makeShape2DList $
  [ (Points (centerC, centerR) minDist (IncompleteCircle thetaDeg0 thetaDeg1 r))
  , (Points (centerC, centerR) minDist (Corner thetaDeg0 thetaDeg1 r))
  ]
makeShape2D (Points (centerC, centerR) minDist (KanizsaTriangle1 thetaDeg len0 len1 r)) =
  let up =
        getShape2DIndexList . makeShape2D . Points (0, 0) len0 $
        (ETriangle (30 + thetaDeg) len0)
      down =
        getShape2DIndexList . makeShape2D . Points (0, 0) len0 $
        (ETriangle (210 + thetaDeg) len0)
   in removeDuplicate minDist .
      L.map (\(x, y) -> (x + centerC, y + centerR)) . makeShape2DList $
      (L.map
         (\((x, y), tDeg) ->
            (Points
               (fromIntegral x, fromIntegral y)
               minDist
               (Corner tDeg 60 len1))) .
       L.zip up . L.map (+ thetaDeg) $
       [210, 90, 330]) L.++
      (L.map
         (\((x, y), tDeg) ->
            (Points (fromIntegral x, fromIntegral y) minDist (PacMan tDeg 60 r))) .
       L.zip down . L.map (+ thetaDeg) $
       [30, 270, 150])
makeShape2D (Points (centerC, centerR) minDist (TJunction thetaDeg len)) =
  removeDuplicate minDist .
  makeShape2DList .
  L.map (\t -> Points (centerC, centerR) minDist (Line (t + thetaDeg) len 0)) $
  [0, 270, 180]
makeShape2D (Points (centerC, centerR) minDist (Cross thetaDeg len)) =
  removeDuplicate minDist .
  makeShape2DList .
  L.map (\t -> Points (centerC, centerR) minDist (Line (t + thetaDeg) len 0)) $
  [0, 90, 180, 270]
makeShape2D (Points (centerC, centerR) minDist (Rectangle thetaDeg width height)) = undefined
makeShape2D (Points (centerC, centerR) minDist (Ehrenstein num r len)) =
  let deltaThetaDeg = 360 / fromIntegral num
      circleIdx = makeShape2D (Points (centerC, centerR) undefined (Circle num r))
   in removeDuplicate minDist .
      makeShape2DList .
      L.zipWith
        (\(i, j) k -> (Points (i, j) minDist (Line (k * deltaThetaDeg) len 1)))
        circleIdx $
      [0 .. (fromIntegral num - 1)]
makeShape2D (Points (centerC, centerR) _ (Circle num r)) =
  let deltaThetaRad = 2 * pi / fromIntegral num
   in [ ( fromIntegral r * cos (n * deltaThetaRad) + centerC
        , fromIntegral r * sin (n * deltaThetaRad) + centerR)
      | n <- [0 .. (fromIntegral num - 1)]
      ]

{-# INLINE makeShape2DContrast #-}
makeShape2DContrast :: Int -> Int -> ContrastArea Shape2D -> R.Array D DIM2 Int
makeShape2DContrast rows cols (ContrastArea whiteBackgroundFlag (Line ori length width)) =
  let centerC = div cols 2
      centerR = div rows 2
   in R.map
        (\x ->
           if whiteBackgroundFlag
             then 255 - x
             else x) .
      R.backpermute
        (Z :. cols :. rows)
        (\(Z :. i' :. j') ->
           let i = i' - centerC
               j = j' - centerR
               x =
                 (centerC) +
                 (round $
                  fromIntegral i * cos (deg2Rad ori) -
                  fromIntegral j * sin (deg2Rad ori))
               y =
                 (centerR) +
                 (round $
                  fromIntegral i * sin (deg2Rad ori) +
                  fromIntegral j * cos (deg2Rad ori))
            in (Z :. x :. y)) .
      fromFunction (Z :. cols :. rows) $ \(Z :. i :. j) ->
        if abs (i - centerC) <= (div length 2) &&
           abs (j - centerR) <= (div width 2)
          then 255
          else 0
makeShape2DContrast rows cols (ContrastArea whiteBackgroundFlag (PacMan thetaDeg0 thetaDeg1 radius)) =
  R.map
    (\x ->
       if whiteBackgroundFlag
         then 255 - x
         else x) .
  fromFunction (Z :. cols :. rows) $ \(Z :. i' :. j') ->
    let i = i' - div cols 2
        j = j' - div rows 2
        r = sqrt . fromIntegral $ i ^ 2 + j ^ 2
        theta = angleFunctionDeg (fromIntegral i) (fromIntegral j)
        deltaTheta = theta - thetaDeg0
     in if r <= fromIntegral radius
          then if L.any
                    (\y -> abs (deltaTheta - y) < thetaDeg1 / 2)
                    [0, 360, -360]
                 then 0
                 else 255
          else 0
makeShape2DContrast rows cols s =
  error $ "makeShape2DContrast: shape is not implemented.\n" L.++ show s
