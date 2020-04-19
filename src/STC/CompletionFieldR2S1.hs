{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module STC.CompletionFieldR2S1 where

import           Array.UnboxedArray          as AU
import           Control.Monad               as M
import           Control.Monad.Parallel      as MP
import           Control.Parallel.Strategies
import           Data.Array.Repa             as R
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           DFT.Plan
import           Filter.Utils
import           Image.IO
import           Image.Transform
import           System.FilePath
import           System.Random
import           Text.Printf
import           Types
import           Utils.Array
import           Utils.Time
import           FokkerPlanck.FourierSeries
import FokkerPlanck.Analytic
import qualified Numeric.LinearAlgebra as NL
import qualified Numeric.LinearAlgebra.HMatrix as NL
import            STC.Utils
import FokkerPlanck.BrownianMotion (thetaPlus)
import Utils.Diagram

-- Implementation of [Williams and Jacbos, 1997]

{-# INLINE checkN #-}
checkN :: Int -> Int -> Int
checkN maxN n
  | n < 0 = checkN maxN (n + maxN)
  | n >= maxN = checkN maxN (n - maxN)
  | otherwise = n 

{-# INLINE rotateST #-}
rotateST ::
     (R.Source r Double)
  => R.Array r DIM3 Double
  -> Int
  -> R.Array U DIM3 Double
rotateST arr 0 = computeUnboxedS . delay $ arr
rotateST arr n' =
  let (Z :. nf :. nx :. ny) = extent arr
      n = checkN nf n'
      deg = (-2) * pi / (fromIntegral nf) * fromIntegral n
  in rotate25D deg (fromIntegral $ center nx, fromIntegral $ center ny) $
     R.backpermute
       (extent arr)
       (\(Z :. k :. i :. j) -> (Z :. (mod (k + nf - n) nf) :. i :. j))
       arr

{-# INLINE timeReversal #-}
timeReversal :: (R.Source r e) => Array r DIM3 e -> Array D DIM3 e
timeReversal arr =
  let (Z :. nf :. _ :. _) = extent arr
      n = div nf 2
  in R.backpermute
       (extent arr)
       (\(Z :. k :. i :. j) -> (Z :. (mod (k + n) nf) :. i :. j))
       arr
       
{-# INLINE rotateSTTR #-}
rotateSTTR ::
     (R.Source r Double)
  => R.Array r DIM3 Double
  -> R.Array U DIM3 Double
rotateSTTR arr =
  let (Z :. _ :. nx :. ny) = extent arr
  in rotate25D pi (fromIntegral $ center nx, fromIntegral $ center ny) arr

{-# INLINE makeR2S1Plan #-}
makeR2S1Plan :: (R.Source r e) => DFTPlan -> R.Array r DIM3 e -> IO DFTPlan
makeR2S1Plan oldPlan arr = do
  let (Z :. orientations :. rows :. cols) = extent arr
  vec3D <-
    VS.map (:+ 0) . VS.fromList <$>
    M.replicateM (orientations * rows * cols) randomIO
  lock <- getFFTWLock
  fst <$>
    (dft1dGPlan lock oldPlan [orientations, rows, cols] [1, 2] vec3D >>= \(plan, vec) ->
       idft1dGPlan lock plan [orientations, rows, cols] [1, 2] vec)

{-# INLINE convolve #-}
convolve ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> Array U DIM3 Double
  -> Array r DIM2 (Complex Double)
  -> IO (VU.Vector Double)
convolve dftPlan arr3D arr2DF = do
  let (Z :. orientations :. rows :. cols) = extent arr3D
      planID3D = DFTPlanID DFT1DG [orientations, rows, cols] [1, 2]
      inversePlanID = DFTPlanID IDFT1DG [orientations, rows, cols] [1, 2]
  vec3DF <-
    dftExecute dftPlan planID3D . VU.convert . VU.map (:+ 0) . toUnboxed $ arr3D
  let arr3DF =
        fromUnboxed (Z :. orientations :. rows :. cols) . VS.convert $ vec3DF
  vec <-
    dftExecute
      dftPlan
      inversePlanID
      (VU.convert . toUnboxed . computeS . R.traverse2 arr3DF arr2DF const $ \f3 f2 idx@(Z :. k :. i :. j) ->
         f3 idx * f2 (Z :. i :. j))
  return . VU.map realPart . VS.convert $ vec

{-# INLINE shareWeightST #-}
shareWeightST ::
     (R.Source r Double)
  => DFTPlan
  -> Array U DIM3 Double
  -> Array r DIM3 Double
  -> IO (Array D DIM3 Double)
shareWeightST dftPlan arr arrG = do
  let (Z :. nf :. _ :. _) = extent arrG
      (Z :. orientations :. rows :. cols) = extent arr
  arrF <-
    fmap (fromUnboxed (extent arr) . VS.convert) .
    dftExecute dftPlan (DFTPlanID DFT1DG [orientations, rows, cols] [1, 2]) .
    VU.convert . VU.map (:+ 0) . toUnboxed $
    arr
  xs <-
    MP.mapM
      (\i ->
         convolve
           dftPlan
           (computeS . makeFilter2D . rotateST arrG $ i)
           (R.slice arrF $ (Z :. i :. All :. All)))
      [0 .. nf - 1]
  return . delay . fromUnboxed (extent arrG) . L.foldl1' (VU.zipWith (+)) $ xs

{-# INLINE computeInitialDistribution #-}
computeInitialDistribution ::
     Int -> Int -> Int -> [R2S1RPPoint] -> R.Array U DIM3 Double
computeInitialDistribution xLen yLen numOrientation xs =
  let deltaTheat = 2 * pi / fromIntegral numOrientation
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
        toUnboxedVector .
        AU.accum (+) 0 ((0, xMin, yMin), (numOrientation - 1, xMax, yMax)) .
        L.map
          (\(R2S1RPPoint (x, y, theta, _)) ->
             ( ((floor $ theta * pi / 180 / deltaTheat), x, y)
             , (1 / (fromIntegral . L.length $ xs)))) $
        xs
   in fromUnboxed (Z :. numOrientation :. xLen :. yLen) vec

-- Implementation of [Williams and Zweck, 2003]

{-# INLINE computeBias #-}
computeBias ::
     Int -> Int -> Int -> [R2S1RPPoint] -> R.Array D DIM3 Double
computeBias xLen yLen numOrientation xs =
  let xShift = div xLen 2
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
        toUnboxedVector .
        AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
        L.map (\(R2S1RPPoint (x, y, _, _)) -> ((x, y), 1)) $
        xs
   in R.traverse
        (fromUnboxed (Z :. xLen :. yLen) vec)
        (const (Z :. numOrientation :. xLen :. yLen)) $ \f (Z :. _ :. i :. j) ->
        f (Z :. i :. j)
        

{-# INLINE computeInitialEigenVec #-}
computeInitialEigenVec ::
     Int -> Int -> Int -> [R2S1RPPoint] -> R.Array D DIM3 Double
computeInitialEigenVec xLen yLen numOrientation xs =
  let xShift = div xLen 2
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
        toUnboxedVector .
        AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
        L.map
          (\(R2S1RPPoint (x, y, _, _)) ->
             ((x, y), (1 / (fromIntegral . L.length $ xs)))) $
        xs
   in R.traverse
        (fromUnboxed (Z :. xLen :. yLen) vec)
        (const (Z :. numOrientation :. xLen :. yLen)) $ \f (Z :. _ :. i :. j) ->
        f (Z :. i :. j) / fromIntegral numOrientation
        
{-# INLINE computeInitialEigenVecGaussian #-}
computeInitialEigenVecGaussian ::
     Int -> Int -> Int -> [(Int, Int, Double)] -> R.Array D DIM3 Double
computeInitialEigenVecGaussian xLen yLen numOrientation xs =
  let xShift = div xLen 2
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
        toUnboxedVector .
        AU.accum (+) 0 ((xMin, yMin), (xMax, yMax)) .
        L.map (\(x, y, v) -> ((x, y), v)) $
        xs
  in R.traverse
       (fromUnboxed (Z :. xLen :. yLen) vec)
       (const (Z :. numOrientation :. xLen :. yLen)) $ \f (Z :. _ :. i :. j) ->
       f (Z :. i :. j) 

{-# INLINE eigenVectorR2S1 #-}
eigenVectorR2S1 ::
     (R.Source s Double, R.Source s1 Double, R.Source s2 Double)
  => DFTPlan
  -> FilePath
  -> R.Array s DIM3 Double
  -> R.Array s DIM3 Double
  -> Int
  -> Bool
  -> R.Array s1 DIM3 Double
  -> R.Array s2 DIM3 Double
  -> IO (R.Array D DIM3 Double)
eigenVectorR2S1 plan folderPath filterSource filterSink n writeFlag bias inputSource = do
  let sourceDist = R.zipWith (*) bias inputSource
      -- s = R.sumAllS sourceDist
      s = L.maximum . R.toList $ sourceDist
      sourceVec = computeS $ R.map (/ s) sourceDist
      (Z :. numOrientation :. _ :. _) = extent inputSource
  printCurrentTime (show n)
  sourceEigenVec <- shareWeightST plan sourceVec filterSource
  -- plotThetaDimension folderPath (printf "R2S1_%d_" n) (31,-4) .
  --   R.backpermute
  --     (extent sourceEigenVec)
  --     (\(Z :. k :. i :. j) ->
  --        (Z :.
  --         (mod (k + numOrientation - (div numOrientation 2)) numOrientation :: Int) :.
  --         i :.
  --         j)) $
  --   sourceEigenVec
  -- sinkEigenVec <-
  --   timeReversal <$>
  --   shareWeightST
  --     plan
  --     (computeS . R.zipWith (*) bias $ sourceEigenVec')
  --     filterSink
  -- let sourceEigenVec =
  --       R.zipWith (\x y -> x + 0.2 * y) sourceEigenVec' sinkEigenVec
  when
    writeFlag
    (plotImageRepa
       (folderPath </> printf "Source_%d.png" n)
       (ImageRepa 8 .
        computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map (\(x) -> (abs x) ** (1/3)) . R.sumS . rotate3D $
        sourceEigenVec))
  return sourceEigenVec

powerMethod ::
     (R.Source s1 Double, R.Source s2 Double)
  => DFTPlan
  -> FilePath
  -> (R.Array U DIM3 Double)
  -> Int
  -> Bool
  -> String
  -> Double
  -> (R.Array s1 DIM3 Double)
  -> (R.Array s2 DIM3 Double)
  -> IO (R.Array U DIM3 Double)
powerMethod oldPlan folderPath filterSource numIteration writeFlag idStr threshold bias' initialEigenVec = do
  let (Z :. numOrientation :. xLen :. yLen) = extent filterSource
      filterSink = rotateSTTR filterSource
      bias = delay bias'
  plotImageRepa
    (folderPath </> "Bias.png")
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
     bias)
  plan <- makeR2S1Plan oldPlan filterSource
  source <-
    M.foldM
      (\input n ->
         eigenVectorR2S1
           plan
           folderPath
           filterSource
           filterSink
           n
           writeFlag
           bias
           input)
      (delay initialEigenVec)
      [1 .. numIteration]
  -- print . R.toList . R.slice source $ (Z :. All :. (div xLen 2) :. (div yLen 2))
  -- plotThetaDimension folderPath "R2S1_Source_" (31,-4) .
  --   R.backpermute
  --     (extent source)
  --     (\(Z :. k :. i :. j) ->
  --        (Z :.
  --         (mod (k + numOrientation - (div numOrientation 2)) numOrientation :: Int) :.
  --         i :.
  --         j)) $
  --   source
  sink <- shareWeightST plan (computeS . R.zipWith (*) bias $ source) filterSink
  plotImageRepa
    (folderPath </> "Sink.png")
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
     sink)
  -- plotThetaDimension folderPath "R2S1_Sink_" (31,-4) .
  --   R.backpermute
  --     (extent sink)
  --     (\(Z :. k :. i :. j) ->
  --        (Z :.
  --         (mod (k + numOrientation - (div numOrientation 2)) numOrientation :: Int) :.
  --         i :.
  --         j)) $
  --   sink
  let completion = computeS $ R.zipWith (*) source sink
      completionR2 = R.sumS . rotate3D $ completion
      completionR2Vec = toUnboxed completionR2
      completionMax = VU.maximum completionR2Vec
      completionMin = VU.minimum completionR2Vec
  plotImageRepa
    (folderPath </> printf "Completion%s.png" idStr)
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> sqrt $ (x - completionMin) / (completionMax - completionMin)) $
     completionR2)
  plotThetaDimension folderPath "R2S1_Completion_" (31, -4) .
    R.backpermute
      (extent completion)
      (\(Z :. k :. i :. j) ->
         (Z :.
          (mod (k + numOrientation - (div numOrientation 2)) numOrientation :: Int) :.
          i :.
          j)) $
    completion
  return completion


eigenVectorR2S1Analytic ::
     Int -> Double -> Double -> Double -> Double -> [(Double, Double)] -> IO [[Double]]
eigenVectorR2S1Analytic !oris !sigma !tau !gamma !delta !xs = do
  let !numPoints = L.length xs
      !locArr = fromListUnboxed (Z :. L.length xs) xs
      !deltaTheta = 2 * pi / fromIntegral oris
      !weights =
        [1 / 16, 1 / 8, 1 / 16, 1 / 8, 1 / 4, 1 / 8, 1 / 16, 1 / 8, 1 / 16]
      !origins =
        L.zipWith
          (\w (i, j) -> (i * delta, j * delta, w))
          weights
          [(i, j) | i <- [-1 .. 1], j <- [-1 .. 1]]
      transitionMatrixArr' =
        R.traverse locArr (\_ -> (Z :. numPoints :. oris :. numPoints :. oris)) $ \f (Z :. i1 :. j1 :. i2 :. j2) ->
          let rowIdx = i1 * oris + j1
              colIdx = i2 * oris + j2
          in if rowIdx /= colIdx && i1 /= i2  
             -- then L.foldl'
             --        (\s (i, j, w) ->
             --           s +
             --           w *
             --           computePji
             --             sigma
             --             tau
             --             (R2S1RP
             --                (fst . f $ (Z :. i2))
             --                (snd . f $ (Z :. i2))
             --                (deltaTheta * fromIntegral j2)
             --                gamma)
             --             (R2S1RP
             --                (i + (fst . f $ (Z :. i1)))
             --                (j + (snd . f $ (Z :. i1)))
             --                (deltaTheta * fromIntegral j1)
             --                gamma))
             --        0
             --        origins
             then computePji
                    sigma
                    tau
                    (R2S1RP
                       ((fst . f $ (Z :. i2)))
                       ((snd . f $ (Z :. i2)))
                       (deltaTheta * fromIntegral j2)
                       gamma)
                    (R2S1RP
                       (fst . f $ (Z :. i1))
                       (snd . f $ (Z :. i1))
                       (deltaTheta * fromIntegral j1)
                       gamma)
             else 0
  transitionMatrixArr <- computeUnboxedP transitionMatrixArr'
      -- transitionMatrixArr = computeUnboxedS transitionMatrixArr'
  let !size = numPoints * oris
      transitionMatrix =
        deepSeqArray transitionMatrixArr . (size NL.>< size) . R.toList $
        transitionMatrixArr
      (eigVal, eigVec) = NL.eig transitionMatrix
      pairs =
        L.reverse .
        L.sortOn (realPart . fst) .
        L.filter
          (\(c, _) ->
             let (m, p) = polar c
             in abs (p) < 1e-10) $
        L.zip (NL.toList eigVal) (L.map NL.toList . NL.toColumns $ eigVec)
      eigenVector =
        fromListUnboxed (Z :. numPoints :. oris) .
        L.map magnitude . snd . L.head $
        pairs
  print . L.map (polar . fst) . L.take 10 $ pairs
  print . L.map (polar . fst) . L.take 10 . L.reverse $ pairs
  print . L.maximum . L.map realPart . snd . L.head $ pairs
  print . L.minimum . L.map realPart . snd . L.head $ pairs
  return . L.map (\i -> R.toList . R.slice eigenVector $ (Z :. i :. All)) $
    [0 .. numPoints - 1]
    
computeContourR2S1Analytic ::
     (R.Source s Double)
  => FilePath
  -> (R.Array U DIM3 Double)
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> R.Array s DIM3 Double
  -> [(Double, Double)]
  -> IO (R.Array U DIM3 Double)
computeContourR2S1Analytic !folderPath !filterSource !oris !rows !cols !sigma !tau !gamma !std !delta !bias !xs = do
  eigenVector <-
    eigenVectorR2S1Analytic
      oris
      sigma
      tau
      gamma
      delta
      (L.map (\(a, b) -> (a * delta, b * delta)) xs)
  plotR2S1
    (folderPath </> "EigenSource.eps")
    (fromIntegral rows)
    (fromIntegral cols)
    (fromIntegral rows / 16)
    xs
    eigenVector
  let n = 0
      idx = [-n .. n]
      idx2D =
        L.filter
          (\(i, j) -> i ^ 2 + j ^ 2 <= n ^ 2)
          [(i, j) | i <- idx, j <- idx]
      ys =
        L.concat $
        L.zipWith
          (\(x', y') thetas ->
             L.concatMap
               (\(i, j) ->
                  let x = round $ x' + fromIntegral i :: Int
                      y = round $ y' + fromIntegral j :: Int
                      v =
                        (exp $
                         ((fromIntegral x - x') ^ 2 + (fromIntegral y - y') ^ 2) *
                         delta ^ 2 /
                         (-2 * std ^ 2)) /
                        (std ^ 2 * 2 * pi)
                  in L.zipWith (\k theta -> ((k, x, y), theta)) [0 ..] thetas)
               idx2D)
          xs
          eigenVector
      (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      eigenArr =
        fromUnboxed (Z :. oris :. cols :. rows) .
        toUnboxedVector .
        AU.accum (+) 0 ((0, minC, minR), (oris - 1, maxC, maxR)) $
        ys
      filterSink = rotateSTTR filterSource
  plan <- makeR2S1Plan emptyPlan filterSource
  source <- shareWeightST plan eigenArr filterSource
  plotImageRepa
    (folderPath </> "Source.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 3)) . R.sumS . rotate3D $
     source)
  plotR2S1Array
    (folderPath </> "Source.eps")
    (fromIntegral rows)
    (fromIntegral cols)
    (fromIntegral rows / 16)
    xs
    source
  sink <- shareWeightST plan eigenArr filterSink
  plotImageRepa
    (folderPath </> "Sink.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 3)) . R.sumS . rotate3D $
     sink)
  plotR2S1Array
    (folderPath </> "Sink.eps")
    (fromIntegral rows)
    (fromIntegral cols)
    (fromIntegral rows / 16)
    xs
    sink
  let completion = computeS $ R.zipWith (*) source sink
      completionR2 = R.sumS . rotate3D $ completion
      completionR2Vec = toUnboxed completionR2
      completionMax = VU.maximum completionR2Vec
      completionMin = VU.minimum completionR2Vec
  print (completionMin, completionMax)
  plotImageRepa
    (folderPath </> printf "Completion.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 2)) .
     R.map (\x -> sqrt $ (x - completionMin) / (completionMax - completionMin)) $
     completionR2)
  plotR2S1Array
    (folderPath </> "Completion.eps")
    (fromIntegral rows)
    (fromIntegral cols)
    (fromIntegral rows / 16)
    xs
    completion
  let completion1 = computeUnboxedS $ completion *^ bias
  source1 <- shareWeightST plan completion1 filterSource    
  plotImageRepa
    (folderPath </> "Source1.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 3)) . R.sumS . rotate3D $
     source1)
  sink1 <- shareWeightST plan completion1 filterSink    
  plotImageRepa
    (folderPath </> "Sink1.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 3)) . R.sumS . rotate3D $
     sink1)
  let completion' = source1 *^ sink1
      completionR2' = R.sumS . rotate3D $ completion'
      completionR2Vec' = toUnboxed completionR2'
      completionMax' = VU.maximum completionR2Vec'
      completionMin' = VU.minimum completionR2Vec'
  plotImageRepa
    (folderPath </> printf "Completion1.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 2)) .
     R.map (\x -> sqrt $ (x - completionMin') / (completionMax' - completionMin')) $
     completionR2')
  return completion

computeContourR2S1Tangent ::
     FilePath
  -> (R.Array U DIM3 Double)
  -> Int
  -> Int
  -> Int
  -> Double
  -> [(Double, Double)]
  -> IO (R.Array U DIM3 Double)
computeContourR2S1Tangent !folderPath !filterSource !oris !rows !cols !delta !xs = do
  let (cMin, cMax) = computeRange cols
      (rMin, rMax) = computeRange rows
      !deltaTheta = 2 * pi / fromIntegral oris
      initArr =
        fromUnboxed (Z :. oris :. cols :. rows) .
        toUnboxedVector .
        AU.accum (+) 0 ((0, cMin, rMin), (oris - 1, cMax, rMax)) .
        L.concatMap
          (\(x, y) ->
             let phi = atan2 y x
                 idx1 = floor $ ((phi `thetaPlus` (pi / 2)) + pi) / deltaTheta
                 idx2 = floor $ ((phi `thetaPlus` ((-pi) / 2)) + pi) / deltaTheta
             in [((idx1, round x, round y), 1), ((idx2, round x, round y), 1)]) $
        xs
      filterSink = rotateSTTR filterSource
  plan <- makeR2S1Plan emptyPlan filterSource
  source <- shareWeightST plan initArr filterSource
  plotImageRepa
    (folderPath </> "Source.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 3)) . R.sumS . rotate3D $
     source)
  sink <- shareWeightST plan initArr filterSink
  plotImageRepa
    (folderPath </> "Sink.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 3)) . R.sumS . rotate3D $
     sink)
  let completion = computeS $ R.zipWith (*) source sink
      completionR2 = R.sumS . rotate3D $ completion
      completionR2Vec = toUnboxed completionR2
      completionMax = VU.maximum completionR2Vec
      completionMin = VU.minimum completionR2Vec
  print (completionMin, completionMax)
  plotImageRepa
    (folderPath </> printf "Completion.png")
    (ImageRepa 8 .
     computeS .
     R.extend (Z :. (1 :: Int) :. All :. All) .
     R.map (\x -> (abs x) ** (1 / 2)) .
     R.map (\x -> sqrt $ (x - completionMin) / (completionMax - completionMin)) $
     completionR2)
  return completion
