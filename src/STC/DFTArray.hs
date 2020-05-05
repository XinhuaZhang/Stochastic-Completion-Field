{-# LANGUAGE BangPatterns #-}
module STC.DFTArray where

import           Array.UnboxedArray         as AU
import           Control.Monad              as M
import           Data.Array.Repa            as R
import           Data.Complex
import           Data.List                  as L
import           Data.Vector.Storable       as VS
import           Data.Vector.Unboxed        as VU
import           DFT.Plan
import           FokkerPlanck.FourierSeries
import           Image.IO
import           STC.Utils
import           System.FilePath
import           Text.Printf
import           Utils.Parallel


data DFTArray =
  DFTArray !Int -- rows
           !Int -- cols
           ![Double] -- thetaFreqs
           ![Double] -- rFreqs
           ![VS.Vector (Complex Double)]

{-# INLINE getDFTArrayVector #-}
getDFTArrayVector :: DFTArray -> [VS.Vector (Complex Double)]
getDFTArrayVector (DFTArray _ _ _ _ vecs) = vecs

{-# INLINE dftArrayToRepa #-}
dftArrayToRepa :: DFTArray -> R.Array U DIM4 (Complex Double)
dftArrayToRepa (DFTArray rows cols thetaFreqs rFreqs vecs) =
  fromUnboxed (Z :. (L.length rFreqs) :. (L.length thetaFreqs) :. cols :. rows) .
  VS.convert . VS.concat $vecs

{-# INLINE repaToDFTArray #-}
repaToDFTArray ::
     [Double] -> [Double] -> R.Array U DIM4 (Complex Double) -> DFTArray
repaToDFTArray thetaFreqs rFreqs arr =
  let (Z :. (!numRFreq) :. (!numThetaFreq) :. cols :. rows) = extent arr
  in DFTArray rows cols thetaFreqs rFreqs .
     L.map
       (\(!rf, !tf) ->
          VS.convert . toUnboxed . computeS . R.slice arr $
          (Z :. rf :. tf :. All :. All)) $
     (,) <$> [0 .. numRFreq - 1] <*> [0 .. numThetaFreq - 1]

{-# INLINE parMapDFTArray #-}
parMapDFTArray ::
     (VS.Vector (Complex Double) -> VS.Vector (Complex Double))
  -> DFTArray
  -> DFTArray
parMapDFTArray f (DFTArray rows cols thetaFreqs rFreqs vecs) =
  DFTArray rows cols thetaFreqs rFreqs . parMap rdeepseq f $ vecs
  
{-# INLINE parZipWithDFTArray #-}
parZipWithDFTArray ::
     (VS.Vector (Complex Double) -> a -> VS.Vector (Complex Double))
  -> DFTArray
  -> [a]
  -> DFTArray
parZipWithDFTArray f (DFTArray rows cols thetaFreqs rFreqs vecs) =
  DFTArray rows cols thetaFreqs rFreqs . parZipWith rdeepseq f vecs 

-- {-# INLINE plotDFTArrayPower #-}
plotDFTArrayPower :: FilePath -> Int -> Int -> DFTArray -> IO ()
plotDFTArrayPower !filePath !rows !cols =
  plotImageRepa filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert . -- VS.map sqrt .
  L.foldl1' (VS.zipWith (+)) .
  parMap rdeepseq (VS.map (\x -> (magnitude x) ** (1) )) . getDFTArrayVector

-- {-# INLINE plotDFTArrayThetaR #-}
plotDFTArrayThetaR :: FilePath -> Int -> Int -> [[Complex Double]] -> DFTArray -> IO ()
plotDFTArrayThetaR !filePath !rows !cols !thetaRHarmonics arr = do
  let vecs =
        computeFourierSeriesThetaR thetaRHarmonics .
        getDFTArrayVector $
        arr
      img =
        VS.convert . L.foldl1' (VS.zipWith (+)) . L.map (VS.map (magnitude)) $
        vecs
      m = VU.maximum img
  plotImageRepa filePath .
    ImageRepa 8 .
    -- computeS .
    -- reduceContrast 100 .
    fromUnboxed (Z :. (1 :: Int) :. cols :. rows) -- . VU.map (\x -> (x / m)** 1)
   $
    img
    -- VS.convert . L.foldl1' (VS.zipWith (+)) . L.map (VS.map magnitude) $
    -- vecs
  -- M.zipWithM_
  --   (\i vec ->
  --      plotImageRepa (dropExtension filePath L.++ (printf "_%d.png" i)) .
  --      ImageRepa 8 .
  --      fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  --      VS.convert . VS.map magnitude $
  --      vec)
  --   [0 .. L.length vecs - 1]
  --   vecs


{-# INLINE plotDFTArrayThetaRMag #-}
plotDFTArrayThetaRMag :: FilePath -> Int -> Int -> [[Complex Double]] -> DFTArray -> IO ()
plotDFTArrayThetaRMag !filePath !rows !cols !thetaRHarmonics =
  plotImageRepa filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert .
  -- VS.map sqrt .
  L.foldl1' (VS.zipWith (max)) .
  L.map (VS.map (\x -> (magnitude x)^2)) .
  computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector

{-# INLINE sparseArrayToDFTArray #-}
sparseArrayToDFTArray ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> [(Int, Int)]
  -> [R.Array U DIM2 (Complex Double)]
  -> DFTArray
sparseArrayToDFTArray rows cols thetaFreqs rFreqs xs sparseArray =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      !numThetaFreq = L.length thetaFreqs
      !numRFreq = L.length rFreqs
  in repaToDFTArray thetaFreqs rFreqs .
     fromUnboxed (Z :. numRFreq :. numThetaFreq :. cols :. rows) .
     toUnboxedVector .
     AU.accum
       (+)
       0
       ((0, 0, minC, minR), (numRFreq - 1, numThetaFreq - 1, maxC, maxR)) .
     L.concat .
     L.zipWith
       (\(x, y) arr ->
          R.toList . R.traverse arr id $ \f idx@(Z :. i :. j) ->
            ((i, j, x, y), f idx))
       xs -- .
     -- L.map
     --   (\arr ->
     --      let !s =
     --            VU.maximum . toUnboxed . computeS . R.map (\x -> (magnitude x)) $
     --            arr
     --      in R.map (/ (s :+ 0)) arr)
  $
     sparseArray
     
{-# INLINE sparseArrayToDFTArray' #-}
sparseArrayToDFTArray' ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> [(Double, Double)]
  -> [R.Array U DIM1 (Complex Double)]
  -> DFTArray
sparseArrayToDFTArray' rows cols thetaFreqs rFreqs xs sparseArray =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      !numThetaFreq = L.length thetaFreqs
      arr =
        fromUnboxed (Z :. (1 :: Int) :. numThetaFreq :. cols :. rows) .
        toUnboxedVector .
        AU.accum (+) 0 ((0, minC, minR), (numThetaFreq - 1, maxC, maxR)) .
        L.concat .
        L.zipWith
          (\(x, y) arr ->
             R.toList . R.traverse arr id $ \f idx@(Z :. j) ->
               ((j, round x, round y), f idx))
          xs -- .
        -- L.map
        --   (\arr ->
        --      let !s =
        --            VU.maximum .
        --            toUnboxed . computeS . R.map (\x -> (magnitude x)) $
        --            arr
        --      in R.map (/ (s :+ 0)) arr)
        $
        sparseArray
  in repaToDFTArray thetaFreqs rFreqs arr 
  

{-# INLINE sparseArrayToDFTArrayG' #-}
sparseArrayToDFTArrayG' ::
     Int
  -> Double
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> [(Double, Double)]
  -> [R.Array U DIM1 (Complex Double)]
  -> DFTArray
sparseArrayToDFTArrayG' n std rows cols thetaFreqs rFreqs xs sparseArray =
  let idx = [-n .. n]
      idx2D =
        L.filter
          (\(i, j) -> i ^ 2 + j ^ 2 < n ^ 2)
          [(i, j) | i <- idx, j <- idx]
      ys =
        L.map
          (\(x', y') ->
             L.map
               (\(i, j) ->
                  let x = round $ x' + fromIntegral i :: Int
                      y = round $ y' + fromIntegral j :: Int
                      v =
                        (exp $
                         ((fromIntegral x - x') ^ 2 + (fromIntegral y - y') ^ 2) /
                         (-2 * std ^ 2)) /
                        (std ^ 2 * 2 * pi)
                  in (x, y, v))
               idx2D)
          xs
      (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      !numThetaFreq = L.length thetaFreqs
      arr =
        fromUnboxed (Z :. (1 :: Int) :. numThetaFreq :. cols :. rows) .
        toUnboxedVector .
        AU.accum (+) 0 ((0, minC, minR), (numThetaFreq - 1, maxC, maxR)) .
        L.concat .
        L.zipWith
          (\zs arr ->
             L.concatMap
               (\(x, y, v) ->
                  R.toList . R.traverse arr id $ \f idx@(Z :. j) ->
                    ((j, x, y), (v :+ 0) * f idx))
               zs)
          ys $
        sparseArray
  in repaToDFTArray thetaFreqs rFreqs arr 


{-# INLINE sparseArrayToDFTArray''' #-}
sparseArrayToDFTArray''' ::
     Int
  -> Int
  -> [Double]
  -> [Double]
  -> [(Double, Double)]
  -> [R.Array U DIM1 (Complex Double)]
  -> [(DFTArray,DFTArray)]
sparseArrayToDFTArray''' rows cols thetaFreqs rFreqs xs sparseArray =
  let (!minR, !maxR) = computeRange rows
      (!minC, !maxC) = computeRange cols
      !numThetaFreq = L.length thetaFreqs
      !pairs =
        L.zip xs -- .
        -- L.map
        --   (\arr ->
        --      let !s =
        --            sqrt . VU.sum .
        --            toUnboxed . computeS . R.map (\x -> (magnitude x)^2) $
        --            arr
        --      in R.map (/ (s :+ 0)) arr) $
        sparseArray
      computeArray pair =
        fromUnboxed (Z :. (1 :: Int) :. numThetaFreq :. cols :. rows) .
        toUnboxedVector .
        AU.accum (+) 0 ((0, minC, minR), (numThetaFreq - 1, maxC, maxR)) .
        L.concatMap
          (\((x, y), arr) ->
             R.toList . R.traverse arr id $ \f idx@(Z :. j) ->
               ((j, round x, round y), f idx)) $
        pair
      func _ [] = []
      func zs (y:ys) =
        (computeArray [y], computeArray (zs L.++ ys)) : func (y : zs) ys
  in L.map
       (\(x, y) ->
          ( repaToDFTArray thetaFreqs rFreqs x
          , repaToDFTArray thetaFreqs rFreqs y)) .
     func [] $
     pairs
