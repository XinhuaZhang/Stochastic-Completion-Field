{-# LANGUAGE BangPatterns #-}
module STC.InitialDistribution where

import           Array.UnboxedArray  as AU
import           Data.Complex
import           Data.List           as L
import           Data.Vector.Generic as VG
import           STC.Convolution
import           STC.DFTArray
import           STC.Point
import           STC.Utils
import           Utils.Parallel

{-# INLINE computeInitialDistribution #-}
computeInitialDistribution ::
     Int -> Int -> [Double] -> [Double] -> Double -> [Point] -> DFTArray
computeInitialDistribution rows cols thetaFreqs rFreqs halfLogPeriod points =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
  in if L.any (\(Point _ _ _ s) -> s <= 0) points
       then error $ "computeInitialDistribution: initial scale <= 0."
       else DFTArray rows cols thetaFreqs rFreqs .
            parMap
              rdeepseq
              (\(!rFreq, !thetaFreq) ->
                 VG.convert .
                 toUnboxedVector .
                 AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
                 L.map
                   (\(Point x y theta scale) ->
                      ( (x, y)
                      , cis $
                        (-pi * log scale) * rFreq / halfLogPeriod -
                        (thetaFreq * theta * pi / 180))) $
                 points) $
            (,) <$> rFreqs <*> thetaFreqs


{-# INLINE computeInitialDistributionPowerMethod #-}
computeInitialDistributionPowerMethod ::
     Int -> Int -> [Double] -> [Double] -> Double -> [Point] -> DFTArray
computeInitialDistributionPowerMethod rows cols thetaFreqs rFreqs halfLogPeriod points =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      !zeroVec = VG.replicate (rows * cols) 0
  in if L.any (\(Point _ _ _ s) -> s <= 0) points
       then error $ "computeInitialDistributionPowerMethod: initial scale <= 0."
       else DFTArray rows cols thetaFreqs rFreqs .
            L.map
              (\(!rFreq, !thetaFreq) ->
                 if -- rFreq == 0 &&
                 thetaFreq == 0
                   then VG.convert .
                        toUnboxedVector .
                        AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
                        L.map (\(Point x y theta scale) -> ((x, y), 1)) $
                        points
                   else zeroVec) $
            (,) <$> rFreqs <*> thetaFreqs
