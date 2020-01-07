{-# LANGUAGE BangPatterns #-}
module STC.CompletionField where

import           Control.Monad.Parallel as MP
import           Data.Array.Repa        as R
import           Data.Complex
import           Data.List              as L
import           Data.Vector.Storable   as VS
import           DFT.Plan
import           Filter.Utils
import           STC.Convolution
import           STC.DFTArray
import           STC.Plan
import           Utils.Parallel

{-# INLINE completionField #-}
completionField :: DFTPlan -> DFTArray -> DFTArray -> IO DFTArray
completionField !plan !source@(DFTArray rows cols thetaFreqs rFreqs _) !sink = do
  let !numThetaFreq = L.length thetaFreqs
      !numRFreq = L.length rFreqs
      !sourceFilter =
        VS.convert .
        toUnboxed .
        computeUnboxedS .
        R.backpermute
          (Z :. numRFreq :. numThetaFreq :. cols :. rows)
          (\(Z :. a :. b :. c :. d) ->
             (Z :. (makeFilterHelper numRFreq a) :.
              (makeFilterHelper numThetaFreq b) :.
              c :.
              d)) .
        dftArrayToRepa $
        source
      !dftID = DFTPlanID DFT1DG [numRFreq, numThetaFreq, cols, rows] [0, 1]
  dftSource <- dftExecute plan dftID $ sourceFilter
  dftSink <- dftExecute plan dftID . VS.concat . getDFTArrayVector $ sink
  arr <-
    fmap
      (fromUnboxed (Z :. numRFreq :. numThetaFreq :. cols :. rows) . VS.convert) .
    dftExecute
      plan
      (DFTPlanID IDFT1DG [numRFreq, numThetaFreq, cols, rows] [0, 1]) .
    VS.zipWith (*) dftSource $
    dftSink
  return $! repaToDFTArray thetaFreqs rFreqs $ arr
