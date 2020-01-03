module FokkerPlanck.Histogramming where

import           Array.UnboxedArray          as UA
import           Data.Complex
import           Data.DList                  as DL
import           Data.List                   as L
import           FokkerPlanck.BrownianMotion
import           FokkerPlanck.Histogram

{-# INLINE sampleScale #-}
sampleScale ::
     Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Particle
  -> Complex Double
sampleScale halfLogPeriod deltaLogRho phiFreq rhoFreq rFreq (Particle phi rho theta r) =
  let n = floor $ (log r + halfLogPeriod) / deltaLogRho
  in (deltaLogRho :+ 0) *
     (L.sum .
      L.map
        (\particle ->
           let (Particle newPhi newRho _ newR) = moveParticle particle
           in cis $
              pi * (rhoFreq * log newRho + rFreq * log newR) / halfLogPeriod +
              phiFreq * newPhi) $
      [(Particle phi rho theta (fromIntegral i * deltaLogRho)) | i <- [1 .. n]])
      



-- {-# INLINE countR2S1 #-}
-- countR2S1 ::
--      (Int, Int)
--   -> (Int, Int)
--   -> Int
--   -> [DList ParticleIndex]
--   -> Histogram Int
-- countR2S1 (xMin, xMax) (yMin, yMax) numOrientations xs =
--   let deltaTheta = 2 * pi / (fromIntegral numOrientations)
--       ys =
--         L.map (\(t, x, y) -> ((floor $ t / deltaTheta, round x, round y), 1)) .
--         DL.toList .
--         DL.concat . L.map (DL.map (\(x, y, t, s, t0, s0) -> (t, x, y))) $
--         xs
--       numTrajectories = L.length ys
--       arr =
--         UA.accum (+) 0 ((0, xMin, yMin), (numOrientations - 1, xMax, yMax)) ys
--    in Histogram
--         [(yMax - yMin + 1), (xMax - xMin + 1), numOrientations]
--         numTrajectories
--         (toUnboxedVector arr)

-- {-# INLINE countR2Z2T0S0 #-}
-- countR2Z2T0S0 ::
--      (Int, Int)
--   -> (Int, Int)
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> [DList ParticleIndex]
--   -> Histogram (Complex Double)
-- countR2Z2T0S0 (xMin, xMax) (yMin, yMax) t0Freqs tFreqs s0Freqs sFreqs maxR xs =
--   let logMaxR =
--         if maxR == 1
--           then 1
--           else log maxR
--       ys =
--         DL.toList .
--         DL.concat .
--         L.map
--           (\((t0f, i), (tf, j), (s0f, k), (sf, l)) ->
--              DL.map
--                (\(x, y, t, s, t0, s0) ->
--                   let !v =
--                         exp
--                           (0 :+
--                            (-t0f * t0 + tf * t +
--                             (sf * s - s0f * s0) * 2 * pi / logMaxR))
--                       !x' = round x
--                       !y' = round y
--                   in ((j, l, i, k, x', y'), v)) .
--              DL.concat $
--              xs) $
--         [ (t0f, tf, s0f, sf)
--         | t0f <- (L.zip t0Freqs [1 ..])
--         , tf <- (L.zip tFreqs [1 ..])
--         , s0f <- (L.zip s0Freqs [1 ..])
--         , sf <- (L.zip sFreqs [1 ..])
--         ]
--       numTrajectories =
--         div
--           (L.length ys)
--           ((L.length t0Freqs) * (L.length tFreqs) * (L.length s0Freqs) *
--            (L.length sFreqs))
--       arr =
--         UA.accum
--           (+)
--           0
--           ( (1, 1, 1, 1, xMin, yMin)
--           , ( L.length tFreqs
--             , L.length sFreqs
--             , L.length t0Freqs
--             , L.length s0Freqs
--             , xMax
--             , yMax))
--           ys
--   in Histogram
--        [ (yMax - yMin + 1)
--        , (xMax - xMin + 1)
--        , (L.length s0Freqs)
--        , (L.length t0Freqs)
--        , (L.length sFreqs)
--        , (L.length tFreqs)
--        ]
--        numTrajectories
--        (toUnboxedVector arr)
