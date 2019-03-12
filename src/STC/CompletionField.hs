{-# LANGUAGE FlexibleContexts #-}

module STC.CompletionField where

import           Array.UnboxedArray        as AU
import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Storable      as VS
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           FokkerPlanck.DomainChange (r2z1Tor2s1)
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           Image.Transform           (normalizeValueRange)
import           System.FilePath           ((</>))
import           System.Random
import           Types
import           Utils.Array

makeR2Z1T0Plan :: DFTPlan -> R2Z1T0Array -> IO DFTPlan
makeR2Z1T0Plan oldPlan arr = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp2 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp3 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [numThetaFreqs, numTheta0Freqs, xLen, yLen]
       [2, 3]
       vecTemp1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [numThetaFreqs, numTheta0Freqs, xLen, yLen]
         [2, 3]
         vec >>= \(plan, _) ->
         dft1dGPlan lock plan [numThetaFreqs, xLen, yLen] [0] vecTemp2 >>= \(plan, vec) ->
           idft1dGPlan lock plan [numThetaFreqs, xLen, yLen] [0] vec >>= \(plan, _) ->
             dft1dGPlan lock plan [numTheta0Freqs, xLen, yLen] [1, 2] vecTemp3 >>= \(plan, vec) ->
               idft1dGPlan lock plan [numTheta0Freqs, xLen, yLen] [1, 2] vec)

computeInitialDistributionR2T0 ::
     DFTPlan -> Int -> Int -> [Double] -> [R2S1RPPoint] -> IO R2T0Array
computeInitialDistributionR2T0 plan xLen yLen theta0Freqs xs =
  let numTheta0Freqs = L.length theta0Freqs
      xShift = div xLen 2
      (xMin, xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      yShift = div yLen 2
      (yMin, yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      vec =
        VU.concat .
        L.map
          (\t0f ->
             toUnboxedVector .
             AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
             L.map
               (\(R2S1RPPoint (x, y, theta, _)) ->
                  ((x, y), exp (0 :+ (-1) * t0f * (theta / 180 * pi)))) $
             xs) $
        theta0Freqs
   in fmap (fromUnboxed (Z :. numTheta0Freqs :. xLen :. yLen) . VS.convert) .
      dftExecute plan (DFTPlanID DFT1DG [numTheta0Freqs, xLen, yLen] [1, 2]) .
      VS.convert $
      vec

{-# INLINE makeFilterR2Z1T0 #-}
makeFilterR2Z1T0 :: R2Z1T0Array -> R2Z1T0Array
makeFilterR2Z1T0 arr =
  let (Z :. _ :. _ :. rows :. cols) = extent arr
   in computeS $
      R.backpermute
        (extent arr)
        (\(Z :. k :. l :. i :. j) ->
           let halfRows = div rows 2
               halfCols = div cols 2
               x =
                 if i < halfRows
                   then i + halfRows
                   else i - halfRows
               y =
                 if j < halfCols
                   then j + halfCols
                   else j - halfCols
            in (Z :. k :. l :. x :. y))
        arr

{-# INLINE makeFilterR2Z1 #-}
makeFilterR2Z1 ::
     (R.Source s (Complex Double))
  => R.Array s DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
makeFilterR2Z1 arr =
  let (Z :. numFreqs :. _ :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k) ->
           let halfNumFreqs = div numFreqs 2
               x =
                 if i < halfNumFreqs
                   then i + halfNumFreqs
                   else i - halfNumFreqs
            in (Z :. x :. j :. k))
        arr

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

{-# INLINE convolveR2T0 #-}
convolveR2T0 :: DFTPlan -> R2Z1T0Array -> R2T0Array -> IO R2Z1T0Array
convolveR2T0 plan filterF initialDistributionF = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent filterF
      arrF =
        computeS $
        R.traverse2
          filterF
          initialDistributionF
          const
          (\f1 f2 idx@(Z :. tf :. t0f :. x :. y) ->
             f1 idx * f2 (Z :. t0f :. x :. y))
  fmap (fromUnboxed (extent filterF) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreqs, numTheta0Freqs, xLen, yLen] [2, 3]) .
    VS.convert . toUnboxed $
    arrF

{-# INLINE timeReverseR2Z1 #-}
timeReverseR2Z1 ::
     (R.Source s (Complex Double))
  => [Double]
  -> R.Array s DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
timeReverseR2Z1 freqs arr =
  R.traverse2
    arr
    (fromListUnboxed (Z :. (L.length freqs)) freqs)
    const
    (\f1 f2 idx@(Z :. k :. i :. j) -> f1 idx * (exp (0 :+ f2 (Z :. k) * pi)))

{-# INLINE timeReversalConvolveR2Z1 #-}
timeReversalConvolveR2Z1 ::
     DFTPlan
  -> [Double]
  -> R.Array U DIM3 (Complex Double)
  -> R.Array U DIM3 (Complex Double)
  -> IO (R.Array U DIM3 (Complex Double))
timeReversalConvolveR2Z1 plan thetaFreqs arr1 arr2 = do
  let (Z :. numThetaFreqs :. xLen :. yLen) = extent arr1
  vec1F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VU.convert . toUnboxed $
    arr1
  vec2F <-
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreqs, xLen, yLen] [0]) .
    VU.convert .
    toUnboxed . computeS . makeFilterR2Z1 . timeReverseR2Z1 thetaFreqs $
    arr2
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute plan (DFTPlanID IDFT1DG [numThetaFreqs, xLen, yLen] [0]) $
    VS.zipWith (*) vec1F vec2F

--- R2Z2T0S0 

makeR2Z2T0S0Plan :: DFTPlan -> R2Z2T0S0Array -> IO DFTPlan
makeR2Z2T0S0Plan oldPlan arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM
      (numThetaFreqs * numScaleFreqs * xLen * yLen * numTheta0Freqs *
       numScale0Freqs)
      randomIO :: IO (VS.Vector (Complex Double))
  vecTemp2 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * numScaleFreqs * xLen * yLen) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp3 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (xLen * yLen * numTheta0Freqs * numScale0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [ numThetaFreqs
       , numScaleFreqs
       , numTheta0Freqs
       , numScale0Freqs
       , xLen
       , yLen
       ]
       [4, 5]
       vecTemp1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [ numThetaFreqs
         , numScaleFreqs
         , numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         ]
         [4, 5]
         vec >>= \(plan, _) ->
         dft1dGPlan
           lock
           plan
           [numThetaFreqs, numScaleFreqs, xLen, yLen]
           [0, 1]
           vecTemp2 >>= \(plan, vec) ->
           idft1dGPlan
             lock
             plan
             [numThetaFreqs, numScaleFreqs, xLen, yLen]
             [0, 1]
             vec >>= \(plan, _) ->
             dft1dGPlan
               lock
               plan
               [numTheta0Freqs, numScale0Freqs, xLen, yLen]
               [2, 3]
               vecTemp3 >>= \(plan, vec) ->
               idft1dGPlan
                 lock
                 plan
                 [numTheta0Freqs, numScale0Freqs, xLen, yLen]
                 [2, 3]
                 vec)

computeInitialDistributionR2T0S0 ::
     DFTPlan
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> [R2S1RPPoint]
  -> IO R2T0S0Array
computeInitialDistributionR2T0S0 plan xLen yLen theta0Freqs scale0Freqs xs =
  let numTheta0Freqs = L.length theta0Freqs
      numScale0Freqs = L.length scale0Freqs
      xShift = div xLen 2
      (xMin, xMax) =
        if odd xLen
          then (-xShift, xShift)
          else (-xShift, xShift - 1)
      yShift = div yLen 2
      (yMin, yMax) =
        if odd yLen
          then (-yShift, yShift)
          else (-yShift, yShift - 1)
      maxR = sqrt . fromIntegral $ xShift ^ 2 + yShift ^ 2
      vec =
        VU.concat .
        L.map
          (\(t0f, s0f) ->
             toUnboxedVector .
             AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
             L.map
               (\(R2S1RPPoint (x, y, theta, scale')) ->
                  let scale =
                        if scale' == 0
                          then 0
                          else log scale'
                   in ( (x, y)
                      , exp
                          (0 :+
                           (-1) *
                           (t0f * (theta / 180 * pi) +
                            s0f * 2 * pi * scale / log maxR)))) $
             xs) $
        [(t0f, s0f) | t0f <- theta0Freqs, s0f <- scale0Freqs]
   in fmap
        (fromUnboxed (Z :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) .
         VS.convert) .
      dftExecute
        plan
        (DFTPlanID DFT1DG [numTheta0Freqs, numScale0Freqs, xLen, yLen] [2, 3]) .
      VS.convert $
      vec

{-# INLINE makeFilterR2Z2T0S0 #-}
makeFilterR2Z2T0S0 :: R2Z2T0S0Array -> R2Z2T0S0Array
makeFilterR2Z2T0S0 arr =
  let (Z :. _ :. _ :. _ :. _ :. rows :. cols) = extent arr
   in computeS $
      R.backpermute
        (extent arr)
        (\(Z :. a :. b :. k :. l :. i :. j) ->
           let halfRows = div rows 2
               halfCols = div cols 2
               x =
                 if i < halfRows
                   then i + halfRows
                   else i - halfRows
               y =
                 if j < halfCols
                   then j + halfCols
                   else j - halfCols
            in (Z :. a :. b :. k :. l :. x :. y))
        arr


{-# INLINE makeFilterR2Z2 #-}
makeFilterR2Z2 ::
     (R.Source s (Complex Double))
  => R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
makeFilterR2Z2 arr =
  let (Z :. numThetaFreqs :. numScaleFreqs :. _ :. _) = extent arr
   in R.backpermute
        (extent arr)
        (\(Z :. i :. j :. k :. l) ->
           let halfNumThetaFreqs = div numThetaFreqs 2
               halfNumScaleFreqs = div numScaleFreqs 2
               x =
                 if i < halfNumThetaFreqs
                   then i + halfNumThetaFreqs
                   else i - halfNumThetaFreqs
               y =
                 if j < halfNumScaleFreqs
                   then j + halfNumScaleFreqs
                   else j - halfNumScaleFreqs
            in (Z :. x :. y :. k :. l))
        arr

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

{-# INLINE convolveR2T0S0 #-}
convolveR2T0S0 :: DFTPlan -> R2Z2T0S0Array -> R2T0S0Array -> IO R2Z2T0S0Array
convolveR2T0S0 plan filterF initialDistributionF = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent filterF
      arrF =
        computeS $
        R.traverse2
          filterF
          initialDistributionF
          const
          (\f1 f2 idx@(Z :. tf :. sf :. t0f :. s0f :. x :. y) ->
             f1 idx * f2 (Z :. t0f :. s0f :. x :. y))
  fmap (fromUnboxed (extent filterF) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID
         IDFT1DG
         [ numThetaFreqs
         , numScaleFreqs
         , numTheta0Freqs
         , numScale0Freqs
         , xLen
         , yLen
         ]
         [4, 5]) .
    VS.convert . toUnboxed $
    arrF

{-# INLINE timeReverseR2Z2 #-}
timeReverseR2Z2 ::
     (R.Source s (Complex Double))
  => [Double]
  -> [Double]
  -> R.Array s DIM4 (Complex Double)
  -> R.Array D DIM4 (Complex Double)
timeReverseR2Z2 thetaFreqs scaleFreqs arr =
  R.traverse2
    arr
    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    const
    (\f1 f2 idx@(Z :. k :. l :. _ :. _) ->
       f1 idx * (exp (0 :+ f2 (Z :. k) * pi)))

{-# INLINE timeReversalConvolveR2Z2 #-}
timeReversalConvolveR2Z2 ::
     DFTPlan
  -> [Double]
  -> [Double]
  -> R.Array U DIM4 (Complex Double)
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
timeReversalConvolveR2Z2 plan thetaFreqs scaleFreqs arr1 arr2 = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen) = extent arr1
  vec1F <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
    VU.convert . toUnboxed $
    arr1
  vec2F <-
    dftExecute
      plan
      (DFTPlanID DFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) .
    VU.convert .
    toUnboxed .
    computeS . makeFilterR2Z2 . timeReverseR2Z2 thetaFreqs scaleFreqs $
    arr2
  fmap (fromUnboxed (extent arr1) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numThetaFreqs, numScaleFreqs, xLen, yLen] [0, 1]) $
    VS.zipWith (*) vec1F vec2F

makeImagePlan ::
     DFTPlan
  -> R.Array U DIM3 Double
  -> IO (DFTPlan, R.Array U DIM3 (Complex Double))
makeImagePlan plan arr = do
  let (Z :. channels :. cols :. rows) = extent arr
  lock <- getFFTWLock
  fmap (\(a, b) -> (a, fromUnboxed (extent arr) . VS.convert $ b)) .
    dft1dGPlan lock plan [channels, cols, rows] [1, 2] .
    VU.convert . VU.map (:+ 0) . toUnboxed $
    arr

{-# INLINE normalizeMagnitude3D #-}
normalizeMagnitude3D ::
     (R.Source r (Complex Double))
  => R.Array r DIM3 (Complex Double)
  -> R.Array D DIM3 (Complex Double)
normalizeMagnitude3D arr =
  let norm =
        R.map sqrt . R.sumS . R.map (\x -> (magnitude x) ^ 2) . rotate3D $ arr
   in R.traverse2
        arr
        norm
        const
        (\f1 f2 idx@(Z :. k :. i :. j) ->
           if f2 (Z :. i :. j) == 0
             then 0
             else f1 idx / (f2 (Z :. i :. j) :+ 0))

-- {-# INLINE normalizeMagnitudeR2Z1 #-}
-- normalizeMagnitudeR2Z1 :: R2T0Array -> R.Array D DIM3 (Complex Double)
-- normalizeMagnitudeR2Z1 arr =
--   let s = R.sumAllS . R.map (\x -> (magnitude x) ^ 2) $ arr
--    in R.map (/ ((sqrt s) :+ 0)) arr
   
{-# INLINE normalizeMagnitudeR2Z1 #-}
normalizeMagnitudeR2Z1 :: R2T0Array -> R.Array D DIM3 (Complex Double)
normalizeMagnitudeR2Z1 arr =
  let s =
        VU.maximum .
        toUnboxed . R.sumS . R.map (\x -> (magnitude x) ^ 2) . rotate3D $
        arr
   in R.map (/ ((sqrt s) :+ 0)) arr

eigenVectorR2Z1Source ::
     DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> R2Z1T0Array
  -> Int
  -> R2T0Array
  -> R2T0Array
  -> IO R2T0Array
eigenVectorR2Z1Source plan folderPath numOrientation thetaFreqs filter n bias input = do
  let (Z :. numThetaFreq :. cols :. rows) = extent input
      maxV =
        VU.maximum .
        VU.map magnitude . toUnboxed . r2z1Tor2s1 numOrientation thetaFreqs $
        input
  initialDistF <-
    fmap (fromUnboxed (Z :. numThetaFreq :. cols :. rows) . VS.convert) .
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreq, cols, rows] [1, 2]) .
    VU.convert .
    toUnboxed . computeS . R.zipWith (*) bias . R.map (/ (maxV :+ 0)) $
    input
  -- Source field
  sourceArr <- convolveR2T0 plan filter initialDistF
  sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceArr
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sourceR2Z1
  plotImageRepaComplex
    (folderPath </> "Source" L.++ "_" L.++ show n L.++ ".png") .
    ImageRepa 8 $
    sourceField
  return sourceR2Z1
  

eigenVectorR2Z1Sink ::
     DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> R2Z1T0Array
  -> Int
  -> R2T0Array
  -> R2T0Array
  -> IO R2T0Array
eigenVectorR2Z1Sink plan folderPath numOrientation thetaFreqs filter n bias input = do
  let (Z :. numThetaFreq :. cols :. rows) = extent input
      maxV =
        VU.maximum .
        VU.map magnitude . toUnboxed . r2z1Tor2s1 numOrientation thetaFreqs $
        input
  initialDistF <-
    fmap (fromUnboxed (Z :. numThetaFreq :. cols :. rows) . VS.convert) .
    dftExecute plan (DFTPlanID DFT1DG [numThetaFreq, cols, rows] [1, 2]) .
    VU.convert .
    toUnboxed . computeS . R.zipWith (*) bias . R.map (/ (maxV :+ 0)) $
    input
  -- Sink field
  sinkArr <- convolveR2T0 plan filter initialDistF
  sinkR2Z1 <- R.sumP . rotateR2Z1T0Array $ sinkArr
  sinkField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sinkR2Z1
  plotImageRepaComplex (folderPath </> "Sink" L.++ "_" L.++ show n L.++ ".png") .
    ImageRepa 8 $
    sinkField
  return sinkR2Z1


completionFieldR2Z1 ::
     DFTPlan -> FilePath -> Int -> [Double]  -> R2T0Array -> R2T0Array -> IO R2T0Array
completionFieldR2Z1 plan folderPath numOrientation thetaFreqs source sink = do
  completionFiled <- timeReversalConvolveR2Z1 plan thetaFreqs source sink
  completionFiledR2 <-
    R.sumP . R.map magnitude . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    completionFiled
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    completionFiledR2
  let avg = R.sumAllS completionFiledR2 / (fromIntegral $ rows * cols)
      (Z :. cols :. rows) = extent completionFiledR2
  plotImageRepa (folderPath </> "Completion_normalized.png") .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . reduceContrast $
    completionFiledR2
  return completionFiled


{-# INLINE reduceContrast #-}
reduceContrast ::
     (R.Source s Double, Shape sh)
  => Array s sh Double
  -> Array D sh Double
reduceContrast arr =
  let avg =
        (R.sumAllS arr) /
        (fromIntegral . L.product . listOfShape . extent $ arr)
   in R.map
        (\y ->
           if y >= avg
             then avg
             else y)
        arr
