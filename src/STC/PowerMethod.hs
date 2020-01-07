{-# LANGUAGE BangPatterns #-}
module STC.PowerMethod where

import           Control.Monad        as M
import           Data.Array.IArray    as IA
import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           DFT.Plan
import           Image.IO
import           STC.CompletionField
import           STC.Convolution
import           STC.DFTArray
import           STC.Plan
import           System.FilePath      ((</>))
import           Text.Printf
import           Utils.Parallel
import           Utils.Time

powerMethod ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> DFTArray
  -> IO DFTArray
powerMethod _ _ _ _ _ _ 0 !arr = return arr
powerMethod !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !numStep !input@(DFTArray rows cols _ _ _) = do
  printCurrentTime (show numStep)
  convolvedArr <-
    convolve Source plan coefficients harmonicsArray input
  let !biasedConvolvedArr = parMapDFTArray (VS.zipWith (*) bias) convolvedArr
      !maxMag =
        L.maximum .
        parMap rdeepseq (VS.maximum . VS.map magnitude) . getDFTArrayVector $
        biasedConvolvedArr
      !normalizedBiasedConvolvedArr =
        parMapDFTArray (VS.map (/ (maxMag :+ 0))) biasedConvolvedArr
  when writeFlag .
    plotDFTArrayPower
      (folderPath </> (printf "Source_%03d.png" numStep))
      cols
      rows $
    convolvedArr
  powerMethod
    plan
    folderPath
    writeFlag
    coefficients
    harmonicsArray
    bias
    (numStep - 1)
    normalizedBiasedConvolvedArr

{-# INLINE computeContour #-}
computeContour ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> DFTArray
  -> IO DFTArray
computeContour !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !numIteration !input@(DFTArray rows cols _ _ _) = do
  eigenVector <-
    powerMethod
      plan
      folderPath
      writeFlag
      coefficients
      harmonicsArray
      bias
      numIteration
      input
  source <- convolve Source plan coefficients harmonicsArray eigenVector
  sink <- convolve Sink plan coefficients harmonicsArray eigenVector
  completion <- completionField plan source sink
  plotDFTArrayPower (folderPath </> "Sink.png") rows cols sink
  plotDFTArrayPower (folderPath </> "CompletionPower.png") rows cols completion
  return completion
