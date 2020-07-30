{-# LANGUAGE Strict #-}
module STC.Bias where

import           Array.UnboxedArray       as AU
import           Data.Array.Repa          as R
import           Data.Complex
import           Data.List                as L
import           Data.Vector.Storable     as VS
import           Data.Vector.Unboxed      as VU
import           DFT.Plan
import           Filter.Utils
import           Pinwheel.FourierSeries2D
import           STC.Point
import           STC.Utils
import           Utils.List

{-# INLINE computeBias #-}
computeBias ::
     (Num a, Unbox a, Storable a) => Int -> Int -> [Point] -> VS.Vector a
computeBias rows cols =
  let (minR, maxR) = computeRange rows
      (minC, maxC) = computeRange cols
  in VS.convert .
     toUnboxedVector .
     AU.accum (+) 0 ((minC, minR), (maxC, maxR)) .
     L.map (\(Point x y theta scale) -> ((round x, round y), 1))

computeBiasPinwheelBasis ::
     DFTPlan -> Int -> Double -> Double -> [Point] -> IO (VS.Vector (Complex Double))
computeBiasPinwheelBasis plan numR2Freq sigma period points =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      envelope =
        computeUnboxedS $
        analyticalFourierCoefficients2 numR2Freq 1 0 0 sigma period (period * sqrt 2)
  in dftExecute plan (DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]) .
     VS.convert .
     toUnboxed .
     computeS .
     makeFilter2D .
     R.zipWith (*) envelope .
     fromListUnboxed (Z :. numR2Freq :. numR2Freq) $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasPinwheelBasis1 ::
     DFTPlan -> Int -> Double -> Double -> [Point] -> VU.Vector (Complex Double)
computeBiasPinwheelBasis1 plan numR2Freq sigma period points =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      envelope =
        toUnboxed . computeS $
        analyticalFourierCoefficients2 numR2Freq 1 0 0 sigma period (period * sqrt 2)
  in VU.zipWith (*) envelope . VU.fromList $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasGaussian ::
     DFTPlan
  -> Int
  -> Double
  -> Double
  -> [Point]
  -> IO (VS.Vector (Complex Double))
computeBiasGaussian plan numR2Freq sigma period points =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      freqCenter = div numR2Freq 2
      envelope =
        fromFunction (Z :. numR2Freq :. numR2Freq) $ \(Z :. i :. j) ->
          (exp $
           (fromIntegral $ (i - freqCenter) ^ 2 + (j - freqCenter) ^ 2) *
           (sigma ^ 2) /
           (-2)) *
          (sigma ^ 2) /
          (2 * pi) :+
          0
  in dftExecute plan (DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]) .
     VS.convert .
     toUnboxed .
     computeS .
     makeFilter2D .
     R.zipWith (*) envelope . fromListUnboxed (Z :. numR2Freq :. numR2Freq) $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]

computeBiasGaussian1 ::
     DFTPlan -> Int -> Double -> Double -> [Point] -> VU.Vector (Complex Double)
computeBiasGaussian1 plan numR2Freq sigma period points =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      freqCenter = div numR2Freq 2
      envelope =
        toUnboxed . computeS . fromFunction (Z :. numR2Freq :. numR2Freq) $ \(Z :. i :. j) ->
          (exp $
           (fromIntegral $ (i - freqCenter) ^ 2 + (j - freqCenter) ^ 2) *
           (sigma ^ 2) /
           (-2)) *
          (sigma ^ 2) /
          (2 * pi) :+
          0
  in VU.zipWith (*) envelope . VU.fromList $
     [ L.foldl'
       (\b (Point x y theta _) ->
          b + (cis $ (freqX * x + freqY * y) * (-2) * pi / period))
       0
       points
     | freqY <- r2Freqs
     , freqX <- r2Freqs
     ]
