{-# LANGUAGE BangPatterns #-}
module STC.DFTArray where

import           Data.Array.Repa            as R
import           Data.Complex
import           Data.List                  as L
import           Data.Vector.Storable       as VS
import           DFT.Plan
import           FokkerPlanck.FourierSeries
import           Image.IO
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

{-# INLINE plotDFTArrayPower #-}
plotDFTArrayPower :: FilePath -> Int -> Int -> DFTArray -> IO ()
plotDFTArrayPower !filePath !rows !cols =
  plotImageRepa filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert .
  L.foldl1' (VS.zipWith (+)) .
  parMap rdeepseq (VS.map (\x -> (magnitude x) ^ 2)) . getDFTArrayVector

{-# INLINE plotDFTArrayThetaR #-}
plotDFTArrayThetaR :: FilePath -> Int -> Int -> [[Complex Double]] -> DFTArray -> IO ()
plotDFTArrayThetaR !filePath !rows !cols !thetaRHarmonics =
  plotImageRepaComplex filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert .
  L.foldl1' (VS.zipWith (+)) .
  computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector

{-# INLINE plotDFTArrayThetaRMag #-}
plotDFTArrayThetaRMag :: FilePath -> Int -> Int -> [[Complex Double]] -> DFTArray -> IO ()
plotDFTArrayThetaRMag !filePath !rows !cols !thetaRHarmonics =
  plotImageRepa filePath .
  ImageRepa 8 .
  fromUnboxed (Z :. (1 :: Int) :. cols :. rows) .
  VS.convert . VS.map sqrt .
  L.foldl1' (VS.zipWith (+)) .
  L.map (VS.map (\x -> (magnitude x)^2)) .
  computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector
