{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}

module STC.CompletionField
  ( module STC.CompletionField
  , module STC.Plan
  , module STC.InitialDistribution
  , module STC.Bias
  , module STC.Convolution
  , module STC.Utils
  ) where

import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Storable      as VS
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           FokkerPlanck.DomainChange (r2z1Tor2s1, r2z2Tor2s1rp)
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           Image.Transform           (normalizeValueRange)
import           STC.Bias
import           STC.Convolution
import           STC.InitialDistribution
import           STC.Plan
import           STC.Utils
import           System.FilePath           ((</>))
import           System.Random
import           Text.Printf
import           Types
import           Utils.Array

{-# INLINE dftR2Z1T0 #-}
dftR2Z1T0 :: DFTPlan -> R2Z1T0Array -> IO R2Z1T0Array
dftR2Z1T0 plan arr = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numTheta0Freqs, xLen, yLen] [2, 3]) .
    VS.convert . toUnboxed $
    arr

{-# INLINE computeSinkFromSourceR2Z1T0 #-}
computeSinkFromSourceR2Z1T0 ::
     (R.Source r (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array r DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
computeSinkFromSourceR2Z1T0 thetaFreqs theta0Freqs sourceArr =
  R.traverse3
    sourceArr
    (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
    (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
    (\a _ _ -> a) $ \f1 f2 f3 idx@(Z :. k :. l :. i :. j) ->
    f1 idx * exp (0 :+ (f2 (Z :. k) + f3 (Z :. l)) * pi)

--- R2Z2T0S0 

{-# INLINE dftR2Z2T0S0 #-}
dftR2Z2T0S0 :: DFTPlan -> R2Z2T0S0Array -> IO R2Z2T0S0Array
dftR2Z2T0S0 plan arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent arr
  fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [ numThetaFreqs
         , numScaleFreqs
         , numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         ]
         [4, 5]) .
    VS.convert . toUnboxed $
    arr

{-# INLINE computeSinkFromSourceR2Z2T0S0 #-}
computeSinkFromSourceR2Z2T0S0 ::
     (R.Source r (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array r DIM6 (Complex Double)
  -> R.Array D DIM6 (Complex Double)
computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceArr =
  R.traverse3
    sourceArr
    (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
    (fromListUnboxed (Z :. L.length theta0Freqs) theta0Freqs)
    (\a _ _ -> a) $ \f1 f2 f3 idx@(Z :. k :. _ :. l :. _ :. i :. j) ->
    f1 idx * exp (0 :+ (f2 (Z :. k) + f3 (Z :. l)) * pi) 

completionFieldR2S1 ::
     DFTPlan
  -> FilePath
  -> String
  -> Int
  -> [Double]
  -> R2T0Array
  -> R2T0Array
  -> IO (R.Array D DIM2 Double)
completionFieldR2S1 plan folderPath idStr numOrientation thetaFreqs source sink = do
  let completionFiled =
        R.zipWith
          (*)
          (r2z1Tor2s1 numOrientation thetaFreqs $ source)
          (r2z1Tor2s1 numOrientation thetaFreqs $ sink)
  completionFiledR2 <- R.sumP . rotate3D $ completionFiled
  plotImageRepaComplex (folderPath </> printf "CompletioR2S1%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionR2S1_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map (logBase 10000) . normalizeValueRange (1, 256) . R.map magnitude $
    completionFiledR2
  let vec =
        toUnboxed .
        computeS .
        R.map (logBase 10000) . normalizeValueRange (1, 256) . R.map magnitude $
        completionFiledR2
  print . L.take 100 . L.reverse . L.sort . VU.toList $ vec
  return . R.map magnitude $ completionFiledR2

completionFieldR2Z1 ::
     DFTPlan -> FilePath -> String -> Int -> [Double]  -> R2T0Array -> R2T0Array -> IO (R.Array D DIM2 Double)
completionFieldR2Z1 plan folderPath idStr numOrientation thetaFreqs source sink = do
  completionFiled <- convolveR2Z1 plan thetaFreqs source sink
  completionFiledR2 <-
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $ completionFiled
  plotImageRepaComplex (folderPath </> printf "Completion%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionMagnitude%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
  let avg = R.sumAllS completionFiledR2 / (fromIntegral $ rows * cols)
      (Z :. cols :. rows) = extent completionFiledR2
  plotImageRepa (folderPath </> printf "Completion_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map log . normalizeValueRange (1, 256) . R.map magnitude $
    completionFiledR2
  return . R.map magnitude $ completionFiledR2

completionFieldR2Z2 ::
     DFTPlan
  -> FilePath
  -> String
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> R2T0S0Array
  -> R2T0S0Array
  -> IO (R.Array D DIM2 Double)
completionFieldR2Z2 plan folderPath idStr numOrientation thetaFreqs numScale scaleFreqs source sink = do
  completionFiled <- convolveR2Z2 plan thetaFreqs scaleFreqs source sink
  completionFiledR2 <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     completionFiled) >>=
    R.sumP
  plotImageRepaComplex (folderPath </> printf "Completion%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  plotImageRepa (folderPath </> printf "CompletionMagnitude%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFiledR2
  let avg = R.sumAllS completionFiledR2 / (fromIntegral $ rows * cols)
      (Z :. cols :. rows) = extent completionFiledR2
  plotImageRepa (folderPath </> printf "Completion_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map log . normalizeValueRange (1, 256) . R.map magnitude $
    completionFiledR2
  return . R.map magnitude $ completionFiledR2

{-# INLINE timeReverseR2Z1T0 #-}
timeReverseR2Z1T0 ::
     (R.Source s (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double) 
timeReverseR2Z1T0 thetaFreqs theta0Freqs arr =
  let newArr =
        R.traverse3
          arr
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          (fromListUnboxed (Z :. (L.length theta0Freqs)) theta0Freqs)
          (\a _ _ -> a)
          (\f1 f2 f3 idx@(Z :. k :. l :. _ :. _) ->
             f1 idx * (exp (0 :+ (-f3 (Z :. l)) * pi)))
   in newArr
