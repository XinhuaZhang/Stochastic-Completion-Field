{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
module FokkerPlanck.Analytic
  ( R2S1RP(..)
  , computePji
  , computePjiCorner
  , computePjiCorner'
  ) where

import           Data.Complex
import           Data.List               as L
import           Numeric.GSL.Polynomials
import           Text.Printf


data R2S1RP = R2S1RP
  { xR2S1RP, yR2S1RP, thetaR2S1RP, gammaR2S1RP :: Double
  } deriving (Show)
             
instance Eq R2S1RP where
  (==) (R2S1RP x1 y1 theta1 _) (R2S1RP x2 y2 theta2 _) =
    (x1 == x2) && (y1 == y2)

{-# INLINE computeCoefficients #-}
computeCoefficients :: R2S1RP -> R2S1RP -> (Double, Double, Double)
computeCoefficients (R2S1RP x_i y_i theta_i _) (R2S1RP x_j y_j theta_j gamma) =
  let x = x_j - x_i
      y = y_j - y_i
      a = (2 + cos (theta_j - theta_i)) / 3
      b =
        (x * (cos theta_j + cos theta_i) + y * (sin theta_j + sin theta_i)) /
        gamma
      c = (x ^ 2 + y ^ 2) / (gamma ^ 2)
   in (a, b, c)

{-# INLINE solveCubicEquation #-}
solveCubicEquation :: Double -> Double -> Double -> Double -> [Complex Double]
solveCubicEquation a b c d =
  let x0 = b ^ 2 - 3 * a * c
      x1 = 2 * b ^ 3 - 9 * a * b * c + 27 * a ^ 2 * d
      x2 = x1 ^ 2 - 4 * x0 ^ 3
      y =
        if x2 >= 0
          then ((x1 + sqrt x2) / 2) ** (1 / 3) :+ 0
          else ((x1 :+ sqrt (abs x2)) / (2 :+ 0)) ** (1 / 3)
      z = (-0.5) :+ 0.8660254037844386 -- (sqrt 3) / 2
  in [ ((b :+ 0) + z ^ k * y + (x0 :+ 0) / (z ^ k * y)) / ((-3) * a :+ 0)
     | k <- [0 .. 2]
     ]

{-# INLINE isRealPositive #-}
isRealPositive :: Complex Double -> Bool
isRealPositive (a :+ b) = (a > 0) && (abs b < 1e-10)

{-# INLINE computePjit #-}
computePjit ::
     Double -> Double -> (Double, Double, Double) -> Double -> Double
computePjit sigma tau (a, b, c) t =
  3 * exp ((-6) * (a * (t ^ 2) - b * t + c) / (sigma * (t ^ 3))) *
  exp (-t / tau) /
  sqrt ((pi * sigma) ^ 3 * (t ^ 7) / 2)

computePji :: Double -> Double -> R2S1RP -> R2S1RP -> Double
computePji sigma tau x_i x_j =
  if x_i == x_j
    then 0
    else let coef@(a, b, c) = computeCoefficients x_i x_j
             roots' =
               polySolve [9 * c / sigma, -6 * b / sigma, 3 * a / sigma, -7 / 4]
             roots = L.filter isRealPositive roots'
             (pOpt, tOpt) =
               L.maximumBy (\x y -> compare (fst x) (fst y)) .
               L.map
                 (\root ->
                    let realRoot = realPart root
                     in (computePjit sigma tau coef realRoot, realRoot)) $
               roots
             f =
               sqrt
                 (2 * pi * (tOpt ^ 5) /
                  (12 * (3 * c - b * tOpt) / sigma + 7 * (tOpt ^ 3) / 2))
          in f * pOpt


{-# INLINE findIntersectionPoint #-}
findIntersectionPoint :: Double -> R2S1RP -> Maybe Double
findIntersectionPoint delta (R2S1RP x y theta _)
  | theta == pi || theta == 0 || z < delta = Nothing
  | otherwise = Just z
  where
    z = x - y / tan theta

{-# INLINE rotate #-}
rotate :: (Double,Double) -> Double -> (Double,Double)
rotate (x, y) theta =
  (x * cos theta - y * sin theta, x * sin theta + y * cos theta)

{-# INLINE computePjiCorner #-}
computePjiCorner :: Double -> Double -> Double -> [Double] -> R2S1RP -> Double
computePjiCorner delta sigma tau orientations point@(R2S1RP x y theta gamma) =
  case findIntersectionPoint delta point of
    Nothing -> 0
    Just z ->
      let p1 = computePji sigma tau (R2S1RP 0 0 0 gamma) (R2S1RP z 0 0 gamma)
                      -- (newX, newY) = rotate (x - z, y) (-theta)
                      -- p2 =
                      --   computePji
                      --     sigma
                      --     tau
                      --     (R2S1RP 0 0 0 gamma)
                      --     (R2S1RP newX newY 0 gamma)
          p2 = computePji sigma tau (R2S1RP z 0 theta gamma) point
                      -- p1 = L.foldl' (\s ori -> s + computePji sigma tau
                      --   (R2S1RP 0 0 0 gamma) (R2S1RP z 0 ori gamma)) 0
                      --   orientations p2 = L.foldl' (\s ori -> s +
                      --   computePji sigma tau (R2S1RP z 0 ori gamma)
                      --   point) 0 orientations
       in p1 * p2


{-# INLINE computePjiCorner' #-}
computePjiCorner' :: Double -> Double -> Double -> Double -> [Double] -> R2S1RP -> Double
computePjiCorner' delta sigma tau threshold orientations point@(R2S1RP x y theta gamma) =
  case findIntersectionPoint delta point of
    Nothing -> 0
    Just z ->
      let p1 = computePji sigma tau (R2S1RP 0 0 0 gamma) (R2S1RP z 0 0 gamma)
       in if p1 < threshold
            then 0
            else let p11 =
                       L.foldl'
                         (\s ori ->
                            s +
                            computePji
                              sigma
                              tau
                              (R2S1RP 0 0 0 gamma)
                              (R2S1RP z 0 ori gamma))
                         0
                         orientations
                                               -- p2 = computePji sigma tau (R2S1RP z 0 theta gamma) point
                     p2 =
                       L.foldl'
                         (\s ori ->
                            s +
                            computePji
                              sigma
                              tau
                              (R2S1RP x y (theta + pi) gamma)
                              (R2S1RP z 0 ori gamma)
                                                -- point
                          )
                         0
                         orientations
                                       -- p3 =
                                       --   L.foldl'
                                       --     (\s ori ->
                                       --        s +
                                       --        computePji
                                       --          sigma
                                       --          tau
                                       --          (R2S1RP 0 0 0 gamma)
                                       --          (R2S1RP z 0 ori gamma) *
                                       --        computePji
                                       --          sigma
                                       --          tau
                                       --          (R2S1RP z 0 theta gamma)
                                       --          point)
                                       --     0
                                       --     orientations
                  in p11 * p2
