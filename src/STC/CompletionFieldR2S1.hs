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
        (\(Z :. k :. i :. j) ->
           let x = nf - n
            in if k <= n - 1
                 then (Z :. k + x :. i :. j)
                 else (Z :. k - n :. i :. j))
        arr

{-# INLINE timeReversal #-}
timeReversal :: (R.Source r e) => Array r DIM3 e -> Array D DIM3 e
timeReversal arr =
  let (Z :. nf :. _ :. _) = extent arr
      n = div nf 2
   in R.backpermute
        (extent arr)
        (\(Z :. k :. i :. j) ->
           let x = nf - n
            in if k >= x
                 then (Z :. k - x :. i :. j)
                 else (Z :. k + n :. i :. j))
        arr
        
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
             ( ((floor $ theta / 360 * fromIntegral numOrientation), x, y)
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
      s = R.sumAllS sourceDist
      sourceVec = computeS $ R.map (/ s) sourceDist
  sourceEigenVec <- shareWeightST plan sourceVec filterSource
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
        computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
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
      filterSink = rotateST filterSource $ div numOrientation 2
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
  sink <- shareWeightST plan (computeS . R.zipWith (*) bias $ source) filterSink
  plotImageRepa
    (folderPath </> "Sink.png")
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
     sink)
  let completion = computeS $ R.zipWith (*) source (timeReversal sink)
  plotImageRepa
    (folderPath </> printf "Completion%s.png" idStr)
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
     completion)
  return completion
