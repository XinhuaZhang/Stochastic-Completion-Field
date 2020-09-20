module GaussianPinwheel where

import           Control.Parallel.Strategies
import           Data.Array.Repa                          as R
import           Data.Complex                             as C
import           Data.List                                as L
import           Data.Vector.Storable                     as VS
import           Data.Vector.Unboxed                      as VU
import           STC.Plan
import           Filter.Utils
import qualified FourierPinwheel.Filtering                as FP
import           FourierPinwheel.GaussianEnvelopePinwheel
import           Graphics.Rendering.Chart.Backend.Cairo
import           Graphics.Rendering.Chart.Easy
import           Image.IO
import           Math.Gamma
import           Pinwheel.FourierSeries2D
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.BLAS
import           Utils.Distribution
import           Utils.List

main = do
  args@(numR2FreqsStr:numThetaFreqStr:numRFreqStr:numPointStr:numThetaStr:numRStr:angularFreqStr:radialFreqStr:stdStr:sigmaStr:_) <-
    getArgs
  let numR2Freqs = read numR2FreqsStr :: Int
      numThetaFreq = read numThetaFreqStr :: Int
      numRFreq = read numRFreqStr :: Int
      numPoint = read numPointStr :: Int
      numTheta = read numThetaStr :: Int
      numR = read numRStr :: Int
      angularFreq = read angularFreqStr :: Int
      radialFreq = read radialFreqStr :: Int
      std = read stdStr :: Double
      sigma = read sigmaStr :: Double
      periodR2 = fromIntegral numR2Freqs
      periodEnv = periodR2 ^ 2 / 4
      folderPath = "output/test/GaussianPinwheel"
  createDirectoryIfMissing True folderPath
  --DFT Plan
  initVec <- generateRadomVector (numR2Freqs ^ 2)
  lock <- getFFTWLock
  plan <-
    fst <$>
    (idft1dGPlan lock emptyPlan [numR2Freqs, numR2Freqs] [0, 1] initVec >>= \(plan, vec) ->
       dft1dGPlan lock plan [numR2Freqs, numR2Freqs] [0, 1] vec)
  let a = 1 / (std * sqrt 2)
      pinwheelCoef =
        centerHollowArray numR2Freqs $
        createFrequencyArray
          numR2Freqs
          (gaussianPinwheelFourierCoefficients
             numR2Freqs
             periodR2
             a
             sigma
             angularFreq
             radialFreq
             periodEnv)
  plotImageRepaComplex (folderPath </> "GaussianPinwheelCoef.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    pinwheelCoef
  pinwheel <-
    dftExecute plan (DFTPlanID IDFT1DG [numR2Freqs, numR2Freqs] [0, 1]) .
    VU.convert . toUnboxed . computeS . makeFilter2D $
    pinwheelCoef
  plotImageRepaComplex (folderPath </> "GaussianPinwheel.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    makeFilter2DInverse .
    fromUnboxed (Z :. numR2Freqs :. numR2Freqs) . VS.convert $
    pinwheel
  let pinwheel1 =
        computeS $
        createFrequencyArray
          numR2Freqs
          (\phi rho ->
             ((exp ((-1) * (rho / std) ^ 2 / 2)) :+ 0) *
             (rho :+ 0) **
             (sigma :+ ((-2) * pi * fromIntegral radialFreq / log periodEnv)) *
             cis ((-phi) * fromIntegral angularFreq))
  pinwheelCoef1 <-
    dftExecute plan (DFTPlanID DFT1DG [numR2Freqs, numR2Freqs] [0, 1]) .
    VU.convert . toUnboxed . computeS . makeFilter2D $
    pinwheel1
  plotImageRepaComplex (folderPath </> "GaussianPinwheelCoef1.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    makeFilter2DInverse .
    fromUnboxed (Z :. numR2Freqs :. numR2Freqs) . VS.convert $
    pinwheelCoef1
  plotImageRepaComplex (folderPath </> "GaussianPinwheel1.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    -- makeFilter2DInverse .
    fromUnboxed (Z :. numR2Freqs :. numR2Freqs) . VS.convert . toUnboxed $
    pinwheel1
