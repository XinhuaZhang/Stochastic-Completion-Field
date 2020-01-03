{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
module STC.Convolution where

import           Control.Arrow
import           Control.Monad.Parallel     as MP (bindM2, mapM)
import           Data.Array.IArray          as IA
import           Data.Array.Repa            as R
import           Data.Complex
import           Data.Ix
import           Data.List                  as L
import           Data.Vector.Storable       as VS
import           Data.Vector.Unboxed        as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.FourierSeries
import           Image.IO
import           Types
import           Utils.Parallel

data Field
  = Source
  | Sink
  deriving (Read, Show)

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

{-# INLINE plotDFTArrayPower #-}
plotDFTArrayPower :: FilePath -> Int -> Int -> DFTArray -> IO ()
plotDFTArrayPower !filePath !rows !cols =
  plotImageRepa filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert .
  VS.map sqrt .
  L.foldl1' (VS.zipWith (+)) .
  parMap rdeepseq (VS.map (\x -> (magnitude x) ^ 2)) . getDFTArrayVector
  
{-# INLINE plotDFTArrayThetaR #-}
plotDFTArrayThetaR :: FilePath -> Int -> Int -> [[Complex Double]] -> DFTArray -> IO ()
plotDFTArrayThetaR filePath rows cols thetaRHarmonics =
  plotImageRepaComplex filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert .
  L.foldl1' (VS.zipWith (+)) .
  computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector
  
{-# INLINE dftHarmonicsArray #-}
dftHarmonicsArray ::
     DFTPlan
  -> Int
  -> Double
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> IO (IA.Array (Int, Int) (VS.Vector (Complex Double)))
dftHarmonicsArray !plan !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff = do
  let !arr =
        computeHarmonicsArray
          numRows
          deltaRow
          numCols
          deltaCol
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          halfLogPeriod
          cutoff
  fmap (listArray (bounds arr)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [numCols, numRows] [0, 1]) .
    L.map
      (VS.convert .
       toUnboxed .
       computeUnboxedS . makeFilter2D . fromUnboxed (Z :. numCols :. numRows)) .
    IA.elems $
    arr

{-# INLINE convolve #-}
convolve ::
     Field
  -> DFTPlan
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> DFTArray
  -> IO DFTArray
convolve !field !plan !coefficients !harmonicsArray !arr@(DFTArray rows cols thetaFreqs rFreqs vecs) = do
  dftVecs <- dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) vecs
  let !initVec = VS.replicate (VS.length . L.head $ vecs) 0
      idx = (,) <$> (L.zip [0 ..] rFreqs) <*> (L.zip [0 ..] thetaFreqs)
  fmap (DFTArray rows cols thetaFreqs rFreqs) .
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) .
    parMap
      rdeepseq
      (\((!r, !rFreq), (!theta, !thetaFreq)) ->
         L.foldl'
           (\(!vec) (((!rho, !rhoFreq), (!phi, !phiFreq)), inputVec) ->
              VS.zipWith
                (+)
                vec
                (VS.map
                   (* (case field of
                         Source ->
                           coefficients R.! (Z :. r :. theta :. rho :. phi)
                         Sink ->
                           (coefficients R.! (Z :. r :. theta :. rho :. phi)) *
                           (cis $ (-(thetaFreq + phiFreq)) * pi))) .
                 VS.zipWith (*) inputVec $
                 (getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq rFreq)))
           initVec .
         L.zip idx $
         dftVecs) $
    idx
