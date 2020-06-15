{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveAnyClass   #-}
{-# LANGUAGE DeriveGeneric    #-}
{-# LANGUAGE FlexibleContexts #-}
module FokkerPlanck.FourierSeries
  ( fourierMellin
  , sampleScale
  -- , computeFourierCoefficients
  , computeHarmonicsArray
  , computeHarmonicsArraySparse
  -- , computeHarmonicsArrayGPU
  , getHarmonics
  , computeThetaRHarmonics
  , computeFourierSeriesThetaR
  -- , computeFourierSeriesR2
  , normalizeFreqArr
  , normalizeFreqArr'
  , plotThetaDimension
  , computeFourierSeriesOfLogPolarHarmonicsArray
  , computeRectangularInverseHarmonics
  ) where

import           Array.UnboxedArray                     as UA
import           Control.DeepSeq
import           Data.Array.IArray                      as IA
import           Data.Array.Repa                        as R
import           Data.Binary
import           Data.Complex
import           Data.DList                             as DL
import           Data.List                              as L
import           Data.Vector.Generic                    as VG
import           Data.Vector.Unboxed                    as VU
import           FokkerPlanck.BrownianMotion
import           FokkerPlanck.Histogram
import           GHC.Generics                           (Generic)
import           Text.Printf
import           Utils.Array
import           Utils.Parallel
import           Utils.Time
import           Graphics.Rendering.Chart.Backend.Cairo
import           Graphics.Rendering.Chart.Easy
import           Image.IO
import           System.FilePath
import           Utils.SimpsonRule
import Numeric.LinearAlgebra as NL
import Numeric.LinearAlgebra.Data as NL
import            Utils.Distribution
import Pinwheel.Base



-- {-# INLINE fourierMellin' #-}
-- fourierMellin' :: Int -> Int -> (Double, Double) -> Complex Double
-- fourierMellin' !angularFreq !radialFreq (!x, !y) =
--   if x == 0 && y == 0
--     then 0
--     else exp $ (-1 * log rho) :+ (tf * atan2 y x + rf * (log rho))

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
       if newRho == 0 -- || newR == 0
         then 0
         else (cos (phiFreq * newPhi + thetaFreq * (theta - newPhi)) :+ 0) *
              (cis $
               (-pi) * (rhoFreq * (newRho) + rFreq * (newR - newRho)) /
               halfLogPeriod))

-- {-# INLINE computeFourierCoefficients #-}
-- computeFourierCoefficients ::
--      [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> Double
--   -> [DList Particle]
--   -> Histogram (Complex Double)
-- computeFourierCoefficients !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !deltaLogRho !xs =
--   Histogram
--     [L.length phiFreqs, L.length rhoFreqs, L.length thetaFreqs, L.length rFreqs]
--     (L.length ys) .
--   toUnboxedVector .
--   UA.accum
--     (+)
--     0
--     ( (1, 1, 1, 1)
--     , ( L.length rFreqs
--       , L.length thetaFreqs
--       , L.length rhoFreqs
--       , L.length phiFreqs)) .
--   L.concat .
--   parMap
--     rdeepseq
--     (\(particle@(Particle phi rho theta r)) ->
--        let !n = floor $ (log r + halfLogPeriod) / deltaLogRho
--            !samples =
--              -- L.map moveParticle  $
--              particle :
--              [ -- (Particle
--              --      phi
--              --      rho
--              --      theta
--              --      (exp $ (fromIntegral i * deltaLogRho - halfLogPeriod)))
--              -- | i <- [0 .. n]
--              ]
--        in L.map
--             (\((rFreq, i), (thetaFreq, j), (rhoFreq, k), (phiFreq, l)) ->
--                ( (i, j, k, l)
--                , -- (deltaLogRho :+ 0) *
--                  1/ (16 * pi^2 * halfLogPeriod :+ 0) *
--                  sampleScale
--                    phiFreq
--                    rhoFreq
--                    thetaFreq
--                    rFreq
--                    halfLogPeriod
--                    samples))
--             freqs) $
--   ys
--   where
--     !ys = DL.toList . DL.concat $ xs
--     !freqs =
--       [ (rFreq', thetaFreq', rhoFreq', phiFreq')
--       | rFreq' <- L.zip rFreqs [1 ..]
--       , thetaFreq' <- L.zip thetaFreqs [1 ..]
--       , rhoFreq' <- L.zip rhoFreqs [1 ..]
--       , phiFreq' <- L.zip phiFreqs [1 ..]
--       ]


{-# INLINE normalizeFreqArr #-}
normalizeFreqArr ::
     Double
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
normalizeFreqArr !std !phiFreqs !rhoFreqs arr =
  computeUnboxedS .
  R.traverse3
    arr
    (fromListUnboxed (Z :. L.length phiFreqs) phiFreqs)
    (fromListUnboxed (Z :. L.length rhoFreqs) rhoFreqs)
    (\sh _ _ -> sh) $ \fArr fPhi fRho idx@(Z :. r :. theta :. rho :. phi) ->
    fArr idx *
    ((exp $
      (-1) *
      ((fPhi (Z :. phi)) ^ 2 + (fPhi (Z :. theta)) ^ 2 + (fRho (Z :. rho)) ^ 2 +
       (fRho (Z :. r)) ^ 2) /
      2 /
      (std ^ 2)) :+
     0)

{-# INLINE normalizeFreqArr' #-}
normalizeFreqArr' ::
     Double
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
normalizeFreqArr' !std !phiFreqs !rhoFreqs arr =
  computeUnboxedS .
  R.traverse3
    arr
    (fromListUnboxed (Z :. L.length phiFreqs) phiFreqs)
    (fromListUnboxed (Z :. L.length rhoFreqs) rhoFreqs)
    (\sh _ _ -> sh) $ \fArr fPhi fRho idx@(Z :. _ :. theta :. rho :. phi) ->
    fArr idx *
    ((exp $
      (-1) * ((fPhi (Z :. phi)) ^ 2 + (fPhi (Z :. theta)) ^ 2) / (2 * std ^ 2) -
      ((fRho (Z :. rho)) ^ 2) / (2 * (std) ^ 2)) :+
     0)

-- {-# INLINE computeHarmonicsArray #-}
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
      rangeFunc1 xs ys = [(L.head xs - L.last ys) .. (L.last xs - L.head ys)]
      rangeFunc2 xs ys =
        (round (L.head xs - L.last ys), round (L.last xs - L.head ys))
      (!tfLB, !tfUB) = rangeFunc2 phiFreqs thetaFreqs
      (!rfLB, !rfUB) = rangeFunc2 rhoFreqs rFreqs
      !xs =
        parMap
          rdeepseq
          (\(!rf, !tf) ->
             let !vec =
                   VG.convert .
                   toUnboxed . computeS . fromFunction (Z :. numCols :. numRows) $ \(Z :. c :. r) ->
                     let !x = fromIntegral (c - centerCol) * deltaCol
                         !y = fromIntegral (r - centerRow) * deltaRow
                         !rho = 0 + (sqrt $ x ^ 2 + y ^ 2)
                         !rho2 =
                           fromIntegral $
                           (c - centerCol) ^ 2 + (r - centerRow) ^ 2
                     in if rho2 > cutoff ^ 2 || (rho <= 0) -- || pi * rho < (abs tf) -- || log ((rho + 1) / rho) > (1 / (2 * rf))
                          then 0
                               -- (x :+ y) ** (tf :+ 0) *
                               -- ((x ^ 2 + y ^ 2) :+ 0) **
                               -- (((-tf - 0.5) :+ rf) / 2)
                          else ((rho :+ 0) ** ((-1) :+ rf)) *
                               (cis (tf * atan2 y x))
             in ((round rf, round tf), vec))
          [ (rf, tf)
          | rf <- rangeFunc1 rhoFreqs rFreqs
          , tf <- rangeFunc1 phiFreqs thetaFreqs
          ]
  in IA.array ((rfLB, tfLB), (rfUB, tfUB)) xs


-- {-# INLINE computeHarmonicsArraySparse #-}
computeHarmonicsArraySparse ::
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
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
computeHarmonicsArraySparse !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff =
  let !centerRow = div numRows 2
      !centerCol = div numCols 2
      rangeFunc1 xs ys =
        [(round $ L.head xs - L.last ys) .. (round $ L.last xs - L.head ys)]
      rangeFunc2 xs ys =
        ((round $ (L.head xs - L.last ys)), (round $ (L.last xs - L.head ys)))
      (!tfLB, !tfUB) = rangeFunc2 phiFreqs thetaFreqs
      (!rfLB, !rfUB) = rangeFunc2 rhoFreqs rFreqs
      !logDR = halfLogPeriod / fromIntegral (min centerRow centerCol)
      !xs =
        parMap
          rseq
          (\(!rf, !tf) ->
             let !arr =
                   computeS . fromFunction (Z :. numCols :. numRows) $ \(Z :. c :. r) ->
                     let !x = fromIntegral (c - centerCol) * deltaCol
                         !y = fromIntegral (r - centerRow) * deltaRow
                         !rho = ((sqrt $ x ^ 2 + y ^ 2))
                         !rho2 =
                           fromIntegral $
                           (c - centerCol) ^ 2 + (r - centerRow) ^ 2
                     in if (rho <= 0) || rho2 > cutoff ^ 2
                          then 0
                          else (x :+ y) ** (fromIntegral tf :+ 0) *
                               (((x ^ 2 + y ^ 2) :+ 0) **
                                (((-(fromIntegral tf) - 1) :+ fromIntegral rf) /
                                 2))
             in deepSeqArray arr ((rf, tf), arr))
          [ (rf, tf)
          | rf <- rangeFunc1 rhoFreqs rFreqs
          , tf <- rangeFunc1 phiFreqs thetaFreqs
          ]
  in IA.array ((rfLB, tfLB), (rfUB, tfUB)) xs

-- {-# INLINE computeHarmonicsArrayGPU #-}
-- computeHarmonicsArrayGPU ::
--      (VG.Vector vector (Complex Double), NFData (vector (Complex Double)))
--   => Int
--   -> Double
--   -> Int
--   -> Double
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> Double
--   -> IA.Array (Int, Int) (vector (Complex Double))
-- computeHarmonicsArrayGPU !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff =
--   let !centerRow = div numRows 2
--       !centerCol = div numCols 2
--       rangeFunc1 xs ys =
--         [round (L.head xs - L.last ys) .. round (L.last xs - L.head ys)]
--       rangeFunc2 xs ys =
--         (round (L.head xs - L.last ys), round (L.last xs - L.head ys))
--       (!tfLB, !tfUB) = rangeFunc2 phiFreqs thetaFreqs
--       (!rfLB, !rfUB) = rangeFunc2 rhoFreqs rFreqs
--       !xs =
--         parMap
--           rdeepseq
--           (\(!rf, !tf) ->
--              let !vec =
--                    VG.convert .
--                    toUnboxed . computeS . fromFunction (Z :. numCols :. numRows) $ \(Z :. c :. r) ->
--                      let !x = fromIntegral (c - centerCol) * deltaCol
--                          !y = fromIntegral (r - centerRow) * deltaRow
--                          !r2 = x ^ 2 + y ^ 2
--                      in if (x == 0 && y == 0) || r2 > cutoff ^ 2
--                           then 0
--                           else cis $
--                                fromIntegral tf * atan2 y x +
--                                (pi * (fromIntegral rf * 0.5 * log r2) /
--                                halfLogPeriod - pi)
--              in ((rf, tf), vec))
--           [ (rf, tf)
--           | rf <-  (rangeFunc1 rhoFreqs rFreqs)
--           , tf <- L.reverse (rangeFunc1 phiFreqs thetaFreqs)
--           ]
--   in IA.array ((rfLB, tfLB), (rfUB, tfUB)) xs

{-# INLINE getHarmonics #-}
getHarmonics ::
     (VG.Vector vector e)
  => IA.Array (Int, Int) (vector e)
  -> Double
  -> Double
  -> Double
  -> Double
  -> vector e
getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq rFreq =
  harmonicsArray IA.! (round (rhoFreq - rFreq), round $ (phiFreq - thetaFreq))


-- {-# INLINE computeFourierSeriesR2 #-}
-- computeFourierSeriesR2 ::
--      Int
--   -> Double
--   -> Int
--   -> Double
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> Double
--   -> IA.Array (Int, Int) (VU.Vector (Complex Double))
--   -> R.Array U DIM4 (Complex Double)
--   -> [VU.Vector (Complex Double)]
-- computeFourierSeriesR2 numRows deltaRow numCols deltaCol phiFreqs rhoFreqs thetaFreqs rFreqs halfLogPeriod cutoff harmonicsArray arr =
--   let (Z :. (!numRFreq) :. (!numThetaFreq) :. _ :. _) = extent arr
--       !initVec = VU.replicate (VU.length (harmonicsArray IA.! (0, 0))) 0
--   in parMap
--        rdeepseq
--        (\((!r, !rFreq), (!theta, !thetaFreq)) ->
--           L.foldl'
--             (\(!vec) ((!rho, !rhoFreq), (!phi, !phiFreq)) ->
--                VU.zipWith
--                  (+)
--                  vec
--                  (VU.map
--                     (* (arr R.! (Z :. r :. theta :. rho :. phi)))
--                     (getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq rFreq)))
--             initVec $
--           (,) <$> (L.zip [0 ..] rhoFreqs) <*> (L.zip [0 ..] phiFreqs)) $
--      (,) <$> (L.zip [0 ..] rFreqs) <*> (L.zip [0 ..] thetaFreqs)

-- {-# INLINE computeThetaRHarmonics #-}
computeThetaRHarmonics ::
     Int -> Int -> [Double] -> [Double] -> Double -> [[(Complex Double)]]
computeThetaRHarmonics !numOrientation !numScale !thetaFreqs !rFreqs !halfLogPeriod =
  let !deltaTheta = 2 * pi / fromIntegral numOrientation
      !deltaScale = 2 * halfLogPeriod / fromIntegral numScale
      !numRFreq = L.length rFreqs
      !numThetaFreq = L.length thetaFreqs
  in if numScale == 1
       then parMap
              rdeepseq
              (\o ->
                 L.map
                   (\freq -> cis $ freq * (fromIntegral o * deltaTheta))
                   thetaFreqs) $
            [0 .. numOrientation - 1]
       else parMap
              rdeepseq
              (\(!o, !s) ->
                 R.toList .
                 R.traverse2
                   (fromListUnboxed (Z :. numRFreq) rFreqs)
                   (fromListUnboxed (Z :. numThetaFreq) thetaFreqs)
                   (\_ _ -> (Z :. numRFreq :. numThetaFreq)) $ \fRFreq fThetaFreq idx@(Z :. rFreq :. thetaFreq) ->
                   -- exp $
                   -- (-0.5 * (fromIntegral s * deltaScale - halfLogPeriod)) :+
                   -- ((fRFreq (Z :. rFreq)) *
                   --  (fromIntegral s * deltaScale - halfLogPeriod) +
                   --  fThetaFreq (Z :. thetaFreq) * fromIntegral o * deltaTheta)
                   cis
                     ((fRFreq (Z :. rFreq)) *
                      (pi * (fromIntegral s * deltaScale - halfLogPeriod) /
                       halfLogPeriod) +
                      fThetaFreq (Z :. thetaFreq) * fromIntegral o * deltaTheta)
               ) $
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
          (\vec0 (!vec1, !v) -> VG.zipWith (+) vec0 . VG.map (* v) $ vec1)
          initVec .
        L.zip vecs) $
     harmonics

plotThetaDimension ::
     (R.Source s Double)
  => FilePath
  -> String
  -> (Int, Int)
  -> R.Array s DIM3 Double
  -> IO ()
plotThetaDimension folderPath prefix (x', y') inputArr = do
  let centerCol = div cols 2
      centerRow = div rows 2
      x = x' + centerCol
      y = y' + centerRow
      (Z :. numOrientation :. cols :. rows) = extent inputArr
      centerArr =
        fromFunction (Z :. (1 :: Int) :. cols :. rows) $ \(Z :. _ :. i :. j) ->
          if i == x && j == y
            then 1 :: Double
            else 0
  plotImageRepa (folderPath </> printf "%s(%d,%d)_center.png" prefix x' y') .
    ImageRepa 8 . computeS $
    centerArr
  let ys' = R.toList . R.slice inputArr $ (Z :. R.All :. x :. y)
      !maxY = L.maximum ys'
      -- ys = L.map (/maxY) ys'
      ys = ys'
      deltaTheta = (360 ::Double) / fromIntegral numOrientation
      xs = [fromIntegral x * deltaTheta | x <- [0 .. numOrientation - 1]]
  toFile def (folderPath </> printf "%s(%d,%d)_theta.png" prefix x' y') $ do
    layout_title .= printf "%s(%d,%d)" prefix x' y'
    layout_x_axis . laxis_generate .= scaledAxis def (0, 359)
    layout_y_axis . laxis_generate .= scaledAxis def (0, maxY)
    plot (line "" [L.zip xs ys])


-- The codes below compute Fourier series in R2

computeFourierSeriesOfLogPolarHarmonicsArray ::
     (VG.Vector vector (Complex Double), NFData (vector (Complex Double)))
  => Double
  -> Double
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> IA.Array (Int, Int) (vector (Complex Double))
computeFourierSeriesOfLogPolarHarmonicsArray !radius !delta !r2Freq !phiFreq !rhoFreq !thetaFreq !rFreq !halfLogPeriod =
  let ((!angularFreqLB, !angularFreqUB), angularFreqs) =
        freqRange phiFreq thetaFreq
      ((!radialFreqLB, !radialFreqUB), radialFreqs) =
        freqRange phiFreq thetaFreq
      !numPoints' = round $ (2 * radius + 1) / delta
      !numPoints =
        if odd numPoints'
          then numPoints'
          else numPoints' - 1
      !center = div numPoints 2
      !period = 2 * radius + 1
      !periodConstant = delta * 2 * pi / period
      !numR2Freq = 2 * r2Freq + 1
      !std = fromIntegral $ div r2Freq 2
      !simpsonWeights =
        toUnboxed . computeS $
        (computeWeightArrFromListOfShape [numPoints, numPoints] :: R.Array D DIM2 (Complex Double))
      !weightedRectangularHarmonics =
        parMap
          rdeepseq
          (\(freq1, freq2) ->
             VU.zipWith (*) simpsonWeights .
             toUnboxed . computeS . fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
               cis
                 (-periodConstant *
                   fromIntegral (freq1 * (i - center) + freq2 * (j - center))))
          [ (freq1, freq2)
          | freq1 <- [-r2Freq .. r2Freq]
          , freq2 <- [-r2Freq .. r2Freq]
          ]
      !cartesianGrid =
        toUnboxed . computeS . fromFunction (Z :. numPoints :. numPoints) $ \(Z :. c :. r) ->
          let !x = fromIntegral (c - center) * delta
              !y = fromIntegral (r - center) * delta
          in (x, y)
      !simpsonNorm = (delta / 3) ^ 2 :+ 0
      !xs =
        parMap
          rdeepseq
          (\(!radialFreq, !angularFreq) ->
             let !logPolarHarmonics =
                   VU.map (fourierMellin 0.5 angularFreq radialFreq) cartesianGrid
                 coefficients =
                   fromListUnboxed (Z :. numR2Freq :. numR2Freq) .
                   L.map (VU.sum . VU.zipWith (*) logPolarHarmonics) $
                   weightedRectangularHarmonics
                 !coefficientsGaussian =
                   VG.convert .
                   toUnboxed . computeS . R.traverse coefficients id $ \f idx@(Z :. xFreq :. yFreq) ->
                     simpsonNorm * (f idx) *
                     (gaussian2D
                        (fromIntegral $ xFreq - r2Freq)
                        (fromIntegral $ yFreq - r2Freq)
                        std :+
                      0)
             in ((radialFreq, angularFreq), coefficientsGaussian))
          [ (radialFreq, angularFreq)
          | radialFreq <- radialFreqs
          , angularFreq <- angularFreqs
          ]
  in IA.array ((radialFreqLB, angularFreqLB), (radialFreqUB, angularFreqUB)) xs

computeRectangularInverseHarmonics ::
     (VG.Vector vector (Complex Double), NFData (vector (Complex Double)))
  => Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> [vector (Complex Double)]
computeRectangularInverseHarmonics !numRows !numCols !delta !radius !maxFreq =
  let !period = 2 * radius + 1
      !periodConstant = delta * 2 * pi / period
      !centerRow = div numRows 2
      !centerCol = div numCols 2
  in parMap
       rdeepseq
       (\(col, row) ->
          VG.fromList
            [ cis
              (periodConstant *
               fromIntegral
                 (freq1 * (row - centerRow) + freq2 * (col - centerCol)))
            | freq1 <- [-maxFreq .. maxFreq]
            , freq2 <- [-maxFreq .. maxFreq]
            ])
       [(col, row) | col <- [0 .. numCols - 1], row <- [0 .. numRows - 1]] 

-- Utilities

{-# INLINE freqRange #-}
freqRange :: Int -> Int -> ((Int, Int), [Int])
freqRange n m =
  let !x = n + m
  in ((-x, x), [-x .. x])
