{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}

module STC.CompletionField where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Storable      as VS
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           FokkerPlanck.DomainChange (r2z1Tor2s1, r2z2Tor2s1rp)
import           Graphics.Gnuplot.Simple
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

{-# INLINE normalizeList #-}
normalizeList :: (Ord e, Fractional e) => [e] -> [e]
normalizeList xs = L.map (/ L.maximum xs) xs

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
    
{-# INLINE dftR2Z2T0S0P #-}
dftR2Z2T0S0P :: DFTPlan -> R2Z2T0S0Array -> IO R2Z2T0S0Array
dftR2Z2T0S0P plan arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent arr
  fromUnboxed (extent arr) . VS.convert . VS.concat <$>
    MP.mapM
      (\(tf, sf) ->
         dftExecute
           plan
           (DFTPlanID
              DFT1DG
              [numTheta0Freqs, numScale0Freqs, xLen, yLen]
              [2, 3]) .
         VS.convert . toUnboxed . computeS . R.slice arr $
         (Z :. tf :. sf :. All :. All :. All :. All))
      [ (tf, sf)
      | tf <- [0 .. numThetaFreqs - 1]
      , sf <- [0 .. numScaleFreqs - 1]
      ]
    
{-# INLINE dftR2Z2T0S0Deconv #-}
dftR2Z2T0S0Deconv :: DFTPlan -> R2Z2T0S0Array -> IO R2Z2T0S0Array
dftR2Z2T0S0Deconv plan arr = do
  let (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen :. numThetaFreqs :. numScaleFreqs) =
        extent arr
  fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [ numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         , numThetaFreqs
         , numScaleFreqs
         ]
         [2, 3, 4, 5]) .
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
  -> IO (R.Array U DIM4 (Complex Double))
completionFieldR2Z2 plan folderPath idStr numOrientation thetaFreqs numScale scaleFreqs source sink = do
  -- completionField <- convolveR2Z2 plan source sink
  completionField <- computeP $ source *^ sink
  let (Z :. _ :. _ :. cols :. rows) = extent completionR2S1RP
      completionR2S1RP =
        r2z2Tor2s1rp
          numOrientation
          thetaFreqs
          numScale
          scaleFreqs
          completionField
      orientationSampleDeg =
        [ 360 * fromIntegral i / fromIntegral numOrientation
        | i <- [0 .. numOrientation - 1]
        ]
  -- MP.mapM_
  --   (\i ->
  --      plotPath
  --        [PNG (folderPath </> printf "Orientation_%d.png" i), Title (show i)] .
  --      L.zip orientationSampleDeg .
  --      normalizeList . R.toList . R.map magnitude . R.slice completionR2S1RP $
  --      (Z :. All :. i :. (div cols 2) :. (div rows 2)))
  --   [0 .. numScale - 1]
  completionFieldR2 <-
    (R.sumP . rotate4D . rotate4D $ completionR2S1RP) >>= R.sumP
  plotImageRepaComplex (folderPath </> printf "Completion%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFieldR2
  (R.sumP . rotate4D . rotate4D $ completionField) >>= R.sumP >>=
    plotImageRepaComplex (folderPath </> printf "CompletionFreq%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All)
  (R.sumP . rotate4D . rotate4D . R.map (\x -> (magnitude x) ^ 2) $
   completionField) >>=
    R.sumP >>=
    plotImageRepa (folderPath </> printf "CompletionFreqMag%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map sqrt
  (R.sumP . rotate4D . rotate4D . R.map (\x -> (magnitude x) ^ 2) $
   completionR2S1RP) >>=
    R.sumP >>=
    plotImageRepa (folderPath </> printf "CompletionMagnitude%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map sqrt
  (R.sumP . rotate4D . rotate4D $ completionR2S1RP) >>= R.sumP >>=
    plotImageRepa (folderPath </> printf "CompletionRealPart%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map realPart
  (R.sumP . rotate4D . rotate4D $ completionR2S1RP) >>= R.sumP >>=
    plotImageRepa
      (folderPath </> printf "CompletionRealPartPositive%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map
      (\x ->
         let y = realPart x
          in if y > 0
               then y
               else 0)
  let avg = R.sumAllS completionFieldR2 / (fromIntegral $ rows * cols)
      (Z :. cols :. rows) = extent completionFieldR2
  plotImageRepa (folderPath </> printf "Completion_normalized%s.png" idStr) .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.map log . normalizeValueRange (1, 256) . R.map magnitude $
    completionFieldR2
  return completionField

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



-- For local eigenvector

{-# INLINE dft4D #-}
dft4D :: DFTPlan -> R2Z2Array -> IO R2Z2Array
dft4D plan arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen) = extent arr
  fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [2, 3]) . 
    VS.convert . toUnboxed $
    arr

{-# INLINE computeSinkFromSourceR2Z2 #-}
computeSinkFromSourceR2Z2 ::
     (R.Source r (Complex Double))
  => [Double]
  -> R.Array r DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
computeSinkFromSourceR2Z2 thetaFreqs sourceArr =
  R.traverse2
    sourceArr
    (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
    (\a _ -> a) $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
    f1 idx * exp (0 :+ f2 (Z :. k) * pi)


completionFieldR2 ::
     DFTPlan
  -> FilePath
  -> String
  -> R2Array
  -> R2Array
  -> IO R2Array
completionFieldR2 plan folderPath idStr source sink = do
  let completionFieldR2 =
        computeS $
        R.zipWith
          (\x y ->  x + y) 
          source
          (R.map
             (\x ->
                let (m, p) = polar x
                 in mkPolar m (p + pi))
             sink)
  plotImageRepaComplex (folderPath </> printf "Completion%s.png" idStr) .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFieldR2
  plotImageRepa (folderPath </> printf "CompletionMagnitude%s.png" idStr) .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map magnitude $
    completionFieldR2
  return completionFieldR2


{-# INLINE timeReverse4D #-}
timeReverse4D ::
     (R.Source s (Complex Double))
  => [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double) 
timeReverse4D thetaFreqs arr =
  let newArr =
        R.traverse2
          arr
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          (\a _ -> a)
          (\f1 f2 idx@(Z :. k :. l :. _ :. _) ->
             f1 idx * (exp (0 :+ (f2 (Z :. k)) * pi)))
   in newArr
   
{-# INLINE flip4D #-}
flip4D ::
     (R.Source s (Complex Double))
  => R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
flip4D arr =
  let (Z :. numThetaFreq :. _ :. _ :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. k :. l :. i :. j) ->
           (Z :. (numThetaFreq - 1 - k) :. l :. i :. j))
        arr
