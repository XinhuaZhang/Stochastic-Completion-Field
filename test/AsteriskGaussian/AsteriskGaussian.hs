{-# LANGUAGE ViewPatterns #-}
module AsteriskGaussian where

import           Control.Parallel.Strategies
import           Data.Array.Repa                        as R
import           Data.Complex                           as C
import           Data.List                              as L
import           Data.Vector.Storable                   as VS
import           Graphics.Rendering.Chart.Backend.Cairo
import           Graphics.Rendering.Chart.Easy
import           Math.Gamma
import           Pinwheel.FourierSeries2D
import           Pinwheel.Gaussian
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.BLAS
import           Utils.List

main = do
  args@(numR2FreqStr:numThetaFreqStr:numPointStr:numThetaStr:xStr:yStr:periodStr:std1Str:std2Str:_) <-
    getArgs
  let numR2Freq = read numR2FreqStr :: Int
      numThetaFreq = read numThetaFreqStr :: Int
      numPoint = read numPointStr :: Int
      numTheta = read numThetaStr :: Int
      x = read xStr :: Double
      y = read yStr :: Double
      period = read periodStr :: Double
      std1 = read std1Str :: Double
      std2 = read std2Str :: Double
      folderPath = "output/test/AsteriskGaussian"
  let centerR2Freq = div numR2Freq 2
      centerThetaFreq = div numThetaFreq 2
      a = std1
      b = std2
      asteriskGaussianFreq =
        VS.concat .
        parMap
          rdeepseq
          (\tFreq' ->
             let tFreq = tFreq' - centerThetaFreq
                 arr =
                   R.map
                     (* ((sqrt $ pi / b) *
                         (exp (fromIntegral tFreq ^ 2 / (-4) / b)) :+
                         0)) $
                   analyticalFourierCoefficients1
                     numR2Freq
                     1
                     tFreq
                     0
                     a
                     period
                     (period * sqrt 2)
             in VS.convert . toUnboxed $
                if tFreq == 0
                  then computeS arr
                  else centerHollowArray numR2Freq arr) $
        [0 .. numThetaFreq - 1] 
      deltaTheta = 2 * pi / (Prelude.fromIntegral numTheta)
      harmonicsD =
        fromFunction (Z :. numTheta :. numThetaFreq :. numR2Freq :. numR2Freq) $ \(Z :. theta' :. tFreq' :. xFreq' :. yFreq') ->
          let tFreq = fromIntegral $ tFreq' - centerThetaFreq
              xFreq = fromIntegral $ xFreq' - centerR2Freq
              yFreq = fromIntegral $ yFreq' - centerR2Freq
              theta = fromIntegral theta' * deltaTheta
          in (1 / (2 * pi * period) :+ 0) *
             (cis ((2 * pi) / period * (x * xFreq + y * yFreq) + tFreq * theta))
  harmonics <- computeUnboxedP harmonicsD
  xs <-
    VS.toList . VS.map magnitude <$>
    gemmBLAS
      numTheta
      1
      (numThetaFreq * numR2Freq ^ 2)
      (VS.convert . toUnboxed $ harmonics)
      asteriskGaussianFreq
  createDirectoryIfMissing True folderPath
  toFile def (folderPath </> "theta.png") $ do
    layout_title .=
      printf "(%d , %d)" (Prelude.round x :: Int) (Prelude.round y :: Int)
    plot
      (line
         ""
         [ L.zip
             [ 2 * (Prelude.fromIntegral i) / (Prelude.fromIntegral numTheta) :: Double
             | i <- [0 .. numTheta - 1]
             ]
             xs
         ])
