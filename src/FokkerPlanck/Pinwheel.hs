module FokkerPlanck.Pinwheel where

import           Data.Array.Repa   as R
import           Data.Complex
import           Data.List         as L
import           Types
import           Utils.Coordinates

{-# INLINE pinwheel #-}
pinwheel :: Double -> Double -> Double -> Double -> Int -> Int -> Complex Double
pinwheel maxR rf af alpha x y
  | r == 0 = 0
  -- | r <= (2 / pi * (abs rf)) = 0
  | otherwise = (((r) :+ 0) ** (alpha :+ (rf * 2 * pi / (log maxR)))) * exp (0 :+ ((af) * theta))
  where
    r = sqrt . fromIntegral $ x ^ (2 :: Int) + y ^ (2 :: Int)
    theta = angleFunctionRad (fromIntegral x) (fromIntegral y)

{-# INLINE computeR2Z2T0S0Array #-}
computeR2Z2T0S0Array ::
     Int
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> IO R2Z2T0S0Array
computeR2Z2T0S0Array xLen yLen alpha thetaFreqs scaleFreqs theta0Freqs scale0Freqs = do
  computeP $
    fromFunction
      (Z :. (L.length thetaFreqs) :. (L.length scaleFreqs) :.
       (L.length theta0Freqs) :.
       (L.length scale0Freqs) :.
       xLen :.
       yLen)
      (\(Z :. tf' :. sf' :. t0f' :. s0f' :. x :. y) ->
         pinwheel
           (sqrt . fromIntegral $ (div xLen 2) ^ 2 + (div yLen 2) ^ 2)
           (fromIntegral (sf' + s0f') + sf + s0f)
           (fromIntegral (tf' + t0f') + tf + t0f)
           alpha
           (x - div xLen 2)
           (y - div yLen 2))
  where
    tf = L.minimum thetaFreqs
    sf = L.minimum scaleFreqs
    t0f = L.minimum theta0Freqs
    s0f = L.minimum scale0Freqs
    
{-# INLINE computeR2Z1T0Array #-}
computeR2Z1T0Array ::
     Int -> Int -> Double -> [Double] -> [Double] -> IO R2Z1T0Array
computeR2Z1T0Array xLen yLen alpha thetaFreqs theta0Freqs =
  computeP $
  fromFunction
    (Z :. (L.length thetaFreqs) :. (L.length theta0Freqs) :. xLen :. yLen)
    (\(Z :. tf' :. t0f' :. x :. y) ->
       pinwheel
         (sqrt . fromIntegral $ (div xLen 2) ^ 2 + (div yLen 2) ^ 2)
         0
         (fromIntegral (tf' + t0f') + tf + t0f)
         alpha
         (x - div xLen 2)
         (y - div yLen 2))
  where
    tf = L.minimum thetaFreqs
    t0f = L.minimum theta0Freqs
