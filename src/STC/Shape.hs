module STC.Shape where

import           Data.Array.Repa    as R
import           Data.Array.Unboxed as AU
import           Data.Ix
import           Data.List          as L

data Points a =
  Points (Double, Double)
         Int
         a
         deriving (Show)

data Shape2D
  = Line { lineOrientationDeg :: Double
         , lineLength         :: Int }
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
  deriving (Show)

{-# INLINE deg2Rad #-}
deg2Rad :: Double -> Double
deg2Rad x = x / 180 * pi

isDuplicate :: (Double,Double) -> [(Double,Double)] -> Bool
isDuplicate (x1, y1) =
  L.any
    (\(x2, y2) ->
       let r = sqrt $ (x1 - x2) ^ 2 + (y1 - y2) ^ 2
        in r <= 2)

{-# INLINE removeDuplicate #-}
removeDuplicate :: [(Double,Double)] -> [(Double,Double)]
removeDuplicate [] = []
removeDuplicate (x:xs) = if isDuplicate x xs
                            then removeDuplicate xs
                            else x : removeDuplicate xs
                            
{-# INLINE getShape2DIndexList #-}
getShape2DIndexList :: [(Double,Double)] -> [(Int, Int)]
getShape2DIndexList = L.map (\(x, y) -> (round x, round y)) . removeDuplicate

{-# INLINE getShape2DIndexListDouble #-}
getShape2DIndexListDouble :: [(Double,Double)] -> [(Double,Double)]
getShape2DIndexListDouble = removeDuplicate

{-# INLINE getShape2DRepaArray #-}
getShape2DRepaArray ::
     Int -> Int -> [(Double,Double)] -> R.Array U DIM3 Double
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
makeShape2D (Points (centerC, centerR) minDist (Line ori len)) =
  L.map (\(x, y) -> (x + centerC, y + centerR)) $
  [ (r * cos (deg2Rad ori), r * sin (deg2Rad ori))
  | r <- (L.map fromIntegral [0,minDist .. len])
  ]
makeShape2D (Points (centerC, centerR) minDist (Corner thetaDeg0 thetaDeg1 len)) =
  let thetaRad0 = deg2Rad thetaDeg0
      thetaRad1 = deg2Rad thetaDeg1
      pointLens = L.map fromIntegral [0,minDist .. len]
   in L.map (\(x, y) -> (x + centerC, y + centerR)) $
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
   in makeShape2DList
        [ (Points
             (x * cos thetaRad, x * sin thetaRad)
             minDist
             (Corner (180 + thetaDeg) 60 len))
        , (Points
             (x * cos theta, x * sin theta)
             minDist
             (Line (90 + thetaDeg) len))
        ]
makeShape2D (Points (centerC, centerR) minDist (IncompleteCircle thetaDeg0 thetaDeg1 r)) =
  let thetaRad0 = deg2Rad thetaDeg0
      thetaRad1 = deg2Rad thetaDeg1
      deltaTheta = 2 * asin (fromIntegral minDist / 2 / fromIntegral r)
   in L.map
        (\(x, y) -> (x + centerC, y + centerR))
        [ (fromIntegral r * cos theta, fromIntegral r * sin theta)
        | theta <-
            [(thetaRad0 + thetaRad1 / 2),(thetaRad0 + thetaRad1 / 2 + deltaTheta) .. (thetaRad0 -
                                                                                      thetaRad1 /
                                                                                      2 +
                                                                                      2 *
                                                                                      pi)]
        ]
makeShape2D (Points (centerC, centerR) minDist (PacMan thetaDeg0 thetaDeg1 r)) =
  makeShape2DList
    [ (Points
         (centerC, centerR)
         minDist
         (IncompleteCircle thetaDeg0 thetaDeg1 r))
    , (Points (centerC, centerR) minDist (Corner thetaDeg0 thetaDeg1 r))
    ]
makeShape2D (Points (centerC, centerR) minDist (KanizsaTriangle1 thetaDeg len0 len1 r)) =
  let up =
        getShape2DIndexList . makeShape2D . Points (0, 0) len0 $
        (ETriangle (30 + thetaDeg) len0)
      down =
        getShape2DIndexList . makeShape2D . Points (0, 0) len0 $
        (ETriangle (210 + thetaDeg) len0)
   in L.map (\(x, y) -> (x + centerC, y + centerR)) . makeShape2DList $
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
