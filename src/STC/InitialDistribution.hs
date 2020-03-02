{-# LANGUAGE BangPatterns #-}
module STC.InitialDistribution where

import           Array.UnboxedArray  as AU
import           Data.Array.Repa     as R
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
computeInitialDistribution rows cols phiFreqs rhoFreqs halfLogPeriod points =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
  in if L.any (\(Point _ _ _ s) -> s <= 0) points
       then error $ "computeInitialDistribution: initial scale <= 0."
       else DFTArray rows cols phiFreqs rhoFreqs .
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
                        (log scale) * (-rFreq) - (thetaFreq * theta * pi / 180))) $
                 points) $
            (,) <$> rhoFreqs <*> phiFreqs


{-# INLINE computeInitialDistribution' #-}
computeInitialDistribution' ::
     Int -> Int -> [Double] -> [Double] -> Double -> [Point] -> DFTArray
computeInitialDistribution' rows cols phiFreqs rhoFreqs  halfLogPeriod points =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
  in if L.any (\(Point _ _ _ s) -> s <= 0) points
       then error $ "computeInitialDistribution: initial scale <= 0."
       else DFTArray rows cols phiFreqs rhoFreqs .
            parMap
              rdeepseq
              (\(!thetaFreq) ->
                 VG.convert .
                 toUnboxedVector .
                 AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
                 L.map
                   (\(Point x y theta scale) ->
                      ((x, y), cis $ -(thetaFreq * theta * pi / 180))) $
                 points) $
            phiFreqs


-- {-# INLINE computeInitialDistributionPowerMethod #-}
-- computeInitialDistributionPowerMethod ::
--      Int -> Int -> [Double] -> [Double] -> Double -> [Point] -> DFTArray
-- computeInitialDistributionPowerMethod rows cols thetaFreqs rFreqs halfLogPeriod points =
--   let (!minR, !maxR) = computeRange rows
--       (!minC, !maxC) = computeRange cols
--       !zeroVec = VG.replicate (rows * cols) 0
--   in if L.any (\(Point _ _ _ s) -> s <= 0) points
--        then error $ "computeInitialDistributionPowerMethod: initial scale <= 0."
--        else DFTArray rows cols thetaFreqs rFreqs .
--             L.map
--               (\(!rFreq, !thetaFreq) ->
--                  if rFreq == 0 && thetaFreq == 0
--                    then VG.convert .
--                         toUnboxedVector .
--                         AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
--                         L.map (\(Point x y theta scale) -> ((x, y), 1)) $
--                         points
--                    else zeroVec) $
--             (,) <$> rFreqs <*> thetaFreqs

{-# INLINE computeInitialDistributionPowerMethod' #-}
computeInitialDistributionPowerMethod' ::
     Int -> Int -> [Double] -> [Double] -> [Point] -> DFTArray
computeInitialDistributionPowerMethod' rows cols phiFreqs rhoFreqs points =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      !zeroVec = VG.replicate (rows * cols) 0
  in if L.any (\(Point _ _ _ s) -> s <= 0) points
       then error $ "computeInitialDistributionPowerMethod: initial scale <= 0."
       else DFTArray rows cols phiFreqs rhoFreqs .
            L.map
              (\(!thetaFreq) ->
                 if thetaFreq == 0
                   then VG.convert .
                        toUnboxedVector .
                        AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
                        L.map (\(Point x y theta scale) -> ((x, y), 1)) $
                        points
                   else zeroVec) $
            phiFreqs


-- {-# INLINE computeInitialDistributionPowerMethodSparse #-}
-- computeInitialDistributionPowerMethodSparse ::
--      [Double] -> [Double] -> [Point] -> [R.Array U DIM2 (Complex Double)]
-- computeInitialDistributionPowerMethodSparse thetaFreqs rFreqs points =
--   let !numThetaFreq = L.length thetaFreqs
--       !numRFreq = L.length rFreqs
--       !thetaCenter = div numThetaFreq 2
--       !rCenter = div numRFreq 2
--       !initArr =
--         computeUnboxedS . fromFunction (Z :. numRFreq :. numThetaFreq) $ \(Z :. i :. j) ->
--           if i == rCenter && j == thetaCenter
--             then 1
--             else 0
--   in L.replicate (L.length points) initArr
                 
  
{-# INLINE computeInitialDistributionPowerMethodSparse' #-}
computeInitialDistributionPowerMethodSparse' ::
     [Double] -> [Double] -> [Point] -> [R.Array U DIM1 (Complex Double)]
computeInitialDistributionPowerMethodSparse' thetaFreqs rFreqs points =
  let !numThetaFreq = L.length thetaFreqs
      !thetaCenter = div numThetaFreq 2
      !initArr =
        computeUnboxedS . fromFunction (Z :. numThetaFreq) $ \(Z :. j) ->
          if j == thetaCenter
            then 1
            else 0
      -- !initArr =
      --   computeUnboxedS .
      --   R.traverse (fromListUnboxed (Z :. numThetaFreq) thetaFreqs) id $ \f (Z :. thetaFreq) ->
      --     cis $ -(fromIntegral thetaFreq * 36 * pi / 180)
  in L.replicate (L.length points) initArr
