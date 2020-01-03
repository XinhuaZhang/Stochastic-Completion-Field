{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveGeneric    #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DeriveAnyClass #-}
module FokkerPlanck.FourierSeries
  ( computeFourierCoefficients
  , computeHarmonicsArray
  , getHarmonics
  , computeThetaRHarmonics
  , computeFourierSeriesThetaR
  , computeFourierSeriesR2
  ) where

import           Array.UnboxedArray          as UA
import           Control.DeepSeq
import           Data.Array.IArray           as IA
import           Data.Array.Repa             as R
import           Data.Binary
import           Data.Complex
import           Data.DList                  as DL
import           Data.List                   as L
import           Data.Vector.Generic         as VG
import           Data.Vector.Unboxed         as VU
import           FokkerPlanck.BrownianMotion
import           FokkerPlanck.Histogram
import           GHC.Generics                (Generic)
import           Text.Printf
import           Utils.Array
import           Utils.Parallel
import           Utils.Time

{-# INLINE sampleScale #-}
sampleScale ::
     Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Particle]
  -> Complex Double
sampleScale !phiFreq !rhoFreq !thetaFreq !rFreq !halfLogPeriod =
  L.sum .
  L.map
    (\(Particle newPhi newRho theta newR) ->
       if newRho == 0 || newR == 0
         then 0
         else (cos (phiFreq * newPhi + thetaFreq * (theta - newPhi)) :+ 0) *
              (cis $
               (-pi) * (rhoFreq * (log newRho) + rFreq * (log (newR / newRho))) /
               halfLogPeriod))

{-# INLINE computeFourierCoefficients #-}
computeFourierCoefficients ::
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> [DList Particle]
  -> Histogram (Complex Double)
computeFourierCoefficients !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !deltaLogRho !xs =
  Histogram
    [L.length phiFreqs, L.length rhoFreqs, L.length thetaFreqs, L.length rFreqs]
    (L.length ys) .
  toUnboxedVector .
  UA.accum
    (+)
    0
    ( (1, 1, 1, 1)
    , ( L.length rFreqs
      , L.length thetaFreqs
      , L.length rhoFreqs
      , L.length phiFreqs)) .
  L.concatMap
    (\(particle@(Particle phi rho theta r)) ->
       let !n = floor $ (log r + halfLogPeriod) / deltaLogRho
           !samples =
             particle :
             (L.map
                moveParticle
                [ (Particle
                     phi
                     rho
                     theta
                     (exp $ (fromIntegral i * deltaLogRho - halfLogPeriod)))
                | i <- [0 .. n]
                ])
       in L.map
            (\((rFreq, i), (thetaFreq, j), (rhoFreq, k), (phiFreq, l)) ->
               ( (i, j, k, l)
               , (deltaLogRho :+ 0) *
                 sampleScale
                   phiFreq
                   rhoFreq
                   thetaFreq
                   rFreq
                   halfLogPeriod
                   samples))
            freqs) $
  ys
  where
    !ys = DL.toList . DL.concat $ xs
    !freqs =
      [ (rFreq', thetaFreq', rhoFreq', phiFreq')
      | rFreq' <- L.zip rFreqs [1 ..]
      , thetaFreq' <- L.zip thetaFreqs [1 ..]
      , rhoFreq' <- L.zip rhoFreqs [1 ..]
      , phiFreq' <- L.zip phiFreqs [1 ..]
      ]

{-# INLINE computeHarmonicsArray #-}
computeHarmonicsArray ::
     (VG.Vector vector (Complex Double), NFData (vector (Complex Double)))
  => Int
  -> Double
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> IA.Array (Int, Int) (vector (Complex Double))
computeHarmonicsArray !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff =
  let !centerRow = div numRows 2
      !centerCol = div numCols 2
      rangeFunc1 xs ys =
        [round (L.head xs - L.last ys) .. round (L.last xs - L.head ys)]
      rangeFunc2 xs ys =
        (round (L.head xs - L.last ys), round (L.last xs - L.head ys))
      (!tfLB, !tfUB) = rangeFunc2 thetaFreqs phiFreqs
      (!rfLB, !rfUB) = rangeFunc2 rFreqs rhoFreqs
      !xs =
        parMap
          rdeepseq
          (\(!tf, !rf) ->
             let !vec =
                   VG.convert .
                   toUnboxed . computeS . fromFunction (Z :. numCols :. numRows) $ \(Z :. c :. r) ->
                     let !x = fromIntegral (c - centerCol) * deltaCol
                         !y = fromIntegral (r - centerRow) * deltaRow
                         !r2 = x ^ 2 + y ^ 2
                     in if (x == 0 && y == 0) || r2 > cutoff ^ 2
                          then 0
                          else cis $
                               fromIntegral tf * atan2 y x +
                               pi * (fromIntegral rf * 0.5 * log r2) /
                               halfLogPeriod
             in ((tf, rf), vec))
          [ (tf, rf)
          | tf <- rangeFunc1 thetaFreqs phiFreqs
          , rf <- rangeFunc1 rFreqs rhoFreqs
          ]
  in IA.array ((tfLB, rfLB), (tfUB, rfUB)) xs
  
{-# INLINE getHarmonics #-}
getHarmonics ::
     (VG.Vector vector (Complex Double))
  => IA.Array (Int, Int) (vector (Complex Double))
  -> Double
  -> Double
  -> Double
  -> Double
  -> vector (Complex Double)
getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq rFreq =
  harmonicsArray IA.! (round (phiFreq - thetaFreq), round (rhoFreq - rFreq))


{-# INLINE computeFourierSeriesR2 #-}
computeFourierSeriesR2 ::
     Int
  -> Double
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> IA.Array (Int, Int) (VU.Vector (Complex Double))
  -> R.Array U DIM4 (Complex Double)
  -> [VU.Vector (Complex Double)]
computeFourierSeriesR2 numRows deltaRow numCols deltaCol phiFreqs rhoFreqs thetaFreqs rFreqs halfLogPeriod cutoff harmonicsArray arr =
  let (Z :. (!numRFreq) :. (!numThetaFreq) :. _ :. _) = extent arr
      !initVec = VU.replicate (VU.length (harmonicsArray IA.! (0, 0))) 0
  in parMap
       rdeepseq
       (\((!r, !rFreq), (!theta, !thetaFreq)) ->
          L.foldl'
            (\(!vec) ((!rho, !rhoFreq), (!phi, !phiFreq)) ->
               VU.zipWith
                 (+)
                 vec
                 (VU.map
                    (* (arr R.! (Z :. r :. theta :. rho :. phi)))
                    (getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq rFreq)))
            initVec $
          (,) <$> (L.zip [0 ..] rhoFreqs) <*> (L.zip [0 ..] phiFreqs)) $
     (,) <$> (L.zip [0 ..] rFreqs) <*> (L.zip [0 ..] thetaFreqs)

{-# INLINE computeThetaRHarmonics #-}
computeThetaRHarmonics ::
     Int -> Int -> [Double] -> [Double] -> Double -> [[(Complex Double)]]
computeThetaRHarmonics !numOrientation !numScale !thetaFreqs !rFreqs !halfLogPeriod =
  let !deltaTheta = 2 * pi / fromIntegral numOrientation
      !deltaScale = 2 * halfLogPeriod / fromIntegral numScale
      !numRFreq = L.length rFreqs
      !numThetaFreq = L.length thetaFreqs
  in parMap
       rdeepseq
       (\(!o, !s) ->
          R.toList .
          R.traverse2
            (fromListUnboxed (Z :. numRFreq) rFreqs)
            (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
            (\_ _ -> (Z :. numRFreq :. numThetaFreq)) $ \fRFreq fThetaFreq idx@(Z :. rFreq :. thetaFreq) ->
            cis
              (pi * (fRFreq (Z :. rFreq)) *
               (fromIntegral s * deltaScale - halfLogPeriod) /
               halfLogPeriod +
               fThetaFreq (Z :. thetaFreq) * fromIntegral o * deltaTheta)) $
     (,) <$> [0 .. numOrientation - 1] <*> [0 .. numScale - 1]


{-# INLINE computeFourierSeriesThetaR #-}
computeFourierSeriesThetaR ::
     (VG.Vector vector (Complex Double), NFData (vector (Complex Double)))
  => [[(Complex Double)]]
  -> [vector (Complex Double)]
  -> [vector (Complex Double)]
computeFourierSeriesThetaR !harmonics !vecs =
  let !initVec = VG.replicate (VG.length . L.head $ vecs) 0
  in parMap
       rdeepseq
       (L.foldl'
          (\(!vec0) (!vec1, !v) -> VG.zipWith (+) vec0 . VG.map (* v) $ vec1)
          initVec .
        L.zip vecs) $
     harmonics
