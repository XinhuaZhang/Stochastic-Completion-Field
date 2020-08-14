{-# LANGUAGE FlexibleContexts #-}
module FourierPinwheel.Filtering where

import           Data.Array.Repa            as R
import           Data.Complex
import           Data.List                  as L
import           Data.Vector.Generic        as VG
import           Data.Vector.Storable       as VS
import           Data.Vector.Unboxed        as VU
import           DFT.Plan
import           Filter.Utils
import           Numeric.GSL.Special.Bessel
import           Utils.List

{-# INLINE idealLowPassFilter #-}
idealLowPassFilter ::
     (VG.Vector vector (Complex Double))
  => Double
  -> Double
  -> Int
  -> Int
  -> vector (Complex Double)
idealLowPassFilter radius period numRhoFreqs numRFreqs =
  let centerRhoFreq = div numRhoFreqs 2
      centerRFreq = div numRFreqs 2
      sinc x =
        if x == 0
          then 1
          else sin (pi * x) / (pi * x)
   in VG.convert .
      toUnboxed . computeS . fromFunction (Z :. numRFreqs :. numRhoFreqs) $ \(Z :. i' :. j') ->
        let i = i' - centerRFreq
            j = j' - centerRhoFreq
         in -- if i == (j)
            --   then
          (2 * radius) / period *
                   sinc (fromIntegral j  / period * (2 * radius)) :+
                   0
              -- else 0

{-# INLINE applyFilterRho #-}
applyFilterRho ::
     DFTPlan
  -> VU.Vector (Complex Double)
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
applyFilterRho plan filterF arr' = do
  let (Z :. numRFreqs :. numThetaFreqs :. numRhoFreqs :. numPhiFreqs) =
        extent arr'
      arr =
        computeUnboxedS $
        R.backpermute
          (Z :. numRFreqs :. numRhoFreqs :. numThetaFreqs :. numPhiFreqs)
          (\(Z :. r :. rho :. theta :. phi) -> (Z :. r :. theta :. rho :. phi))
          arr'
  arrF <-
    fmap (fromUnboxed (extent arr) . VG.convert) .
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         [numRFreqs, numRhoFreqs, numThetaFreqs, numPhiFreqs]
         [0, 1]) .
    VG.convert . toUnboxed $
    arr
  let convolvedArrF =
        computeUnboxedS .
        R.traverse2
          arrF
          (fromUnboxed (Z :. numRFreqs :. numRhoFreqs) filterF)
          const $ \fArr fFilter idx@(Z :. r :. rho :. theta :. phi) ->
          fArr idx * fFilter (Z :. r :. rho)
  fmap
    (computeS .
     R.backpermute
       (Z :. numRFreqs :. numThetaFreqs :. numRhoFreqs :. numPhiFreqs)
       (\(Z :. r :. theta :. rho :. phi) -> (Z :. r :. rho :. theta :. phi)) .
     fromUnboxed (extent arr) . VG.convert) .
    dftExecute
      plan
      (DFTPlanID
         IDFT1DG
         [numRFreqs, numRhoFreqs, numThetaFreqs, numPhiFreqs]
         [0, 1]) .
    VG.convert . toUnboxed $
    convolvedArrF

hollowCoefficients ::
     DFTPlan
  -> Double
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
hollowCoefficients plan radius period coefficients = do
  let (Z :. numRFreqs :. _ :. numRhoFreqs :. _) = extent coefficients
      lpf = idealLowPassFilter (log radius) (log period) numRhoFreqs numRFreqs
      filter =
        VG.convert .
        toUnboxed .
        computeS . makeFilter2D . fromUnboxed (Z :. numRFreqs :. numRhoFreqs) $
        lpf
  filterF <-
    VG.convert <$>
    dftExecute plan (DFTPlanID DFT1DG [numRFreqs, numRhoFreqs] [0, 1]) filter
  -- computeS . R.zipWith (-) coefficients <$>
  applyFilterRho plan filterF coefficients
