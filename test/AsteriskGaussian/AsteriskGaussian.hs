{-# LANGUAGE ViewPatterns #-}
module AsteriskGaussian where

import           Control.Parallel.Strategies
import           Data.Array.Repa                        as R
import           Data.Complex                           as C
import           Data.List                              as L
import           Data.Vector.Storable                   as VS
import           Data.Vector.Unboxed                    as VU
import           Graphics.Rendering.Chart.Backend.Cairo
import           Graphics.Rendering.Chart.Easy
import           Math.Gamma
import           Pinwheel.FourierSeries2D
-- import           Pinwheel.Gaussian
import qualified FourierPinwheel.Filtering              as FP
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.BLAS
import           Utils.List
import Utils.Distribution

main = do
  args@(numR2FreqStr:numThetaFreqStr:numRFreqStr:numPointStr:numThetaStr:numRStr:xStr:yStr:periodStr:std1Str:std2Str:stdR2Str:_) <-
    getArgs
  let numR2Freq = read numR2FreqStr :: Int
      numThetaFreq = read numThetaFreqStr :: Int
      numRFreq = read numRFreqStr :: Int
      numPoint = read numPointStr :: Int
      numTheta = read numThetaStr :: Int
      numR = read numRStr :: Int
      x = read xStr :: Double
      y = read yStr :: Double
      period = read periodStr :: Double
      std1 = read std1Str :: Double
      std2 = read std2Str :: Double
      stdR2 = read stdR2Str :: Double
      folderPath = "output/test/AsteriskGaussian"
  let centerR2Freq = div numR2Freq 2
      centerThetaFreq = div numThetaFreq 2
      periodEnv = period ^ 2 / 4
      a = std1
      b = std2
      asteriskGaussianFreqTheta =
        VS.concat .
        parMap
          rdeepseq
          (\tFreq' ->
             let tFreq = tFreq' - centerThetaFreq
                 arr =
                   R.map (* (gaussian1DFreq (fromIntegral tFreq) b :+ 0)) $
                   analyticalFourierCoefficients1
                     numR2Freq
                     1
                     tFreq
                     0
                     a
                     period
                     periodEnv
              in VS.convert . toUnboxed . computeS $
                 centerHollowArray numR2Freq arr) $
        [0 .. numThetaFreq - 1]
      deltaTheta = 2 * pi / Prelude.fromIntegral numTheta
      harmonicsThetaD =
        fromFunction (Z :. numTheta :. numThetaFreq :. numR2Freq :. numR2Freq) $ \(Z :. theta' :. tFreq' :. xFreq' :. yFreq') ->
          let tFreq = fromIntegral $ tFreq' - centerThetaFreq
              xFreq = fromIntegral $ xFreq' - centerR2Freq
              yFreq = fromIntegral $ yFreq' - centerR2Freq
              theta = fromIntegral theta' * deltaTheta
           in (1 / (2 * pi) :+ 0) *
              cis (2 * pi / period * (x * xFreq + y * yFreq) + tFreq * theta)
  harmonicsTheta <- computeUnboxedP harmonicsThetaD
  xs <-
    VS.toList . VS.map magnitude <$>
    gemmBLAS
      numTheta
      1
      (numThetaFreq * numR2Freq ^ 2)
      (VS.convert . toUnboxed $ harmonicsTheta)
      asteriskGaussianFreqTheta
  toFile def (folderPath </> "theta.png") $ do
    layout_title .=
      printf "(%d , %d)" (Prelude.round x :: Int) (Prelude.round y :: Int)
    plot
      (line
         ""
         [ L.zip
             [fromIntegral i * deltaTheta / pi | i <- [0 .. numTheta - 1]]
             xs
         ])
  let centerRFreq = div numRFreq 2
      centerR2Freq = div numR2Freq 2
      gaussian2D =
        computeUnboxedS . fromFunction (Z :. numR2Freq :. numR2Freq) $ \(Z :. i' :. j') ->
          let i = i' - centerR2Freq
              j = j' - centerR2Freq
           in exp
                (pi * fromIntegral (i ^ 2 + j ^ 2) /
                 ((-1) * period ^ 2 * stdR2 ^ 2)) /
              (2 * pi * stdR2 ^ 2) :+
              0
      asteriskGaussianFreqR =
        VS.concat .
        parMap
          rdeepseq
          (\rFreq' ->
             let rFreq = rFreq' - centerRFreq
                 arr =
                   R.map
                     (* ((gaussian1DFourierCoefficients
                            (fromIntegral rFreq)
                            (log periodEnv)
                            b :+
                          0) -- *
                         -- cis
                         --   ((-2) * pi / log periodEnv * fromIntegral rFreq *
                         --    log 0.25)
                        )) . centerHollowArray numR2Freq $
                   analyticalFourierCoefficients1
                     numR2Freq
                     1
                     0
                     rFreq
                     a
                     period
                     periodEnv
              in VS.convert . toUnboxed . computeS . centerHollowArray numR2Freq $
                 arr -- *^ gaussian2D
           ) $
        [0 .. numRFreq - 1]
      deltaR = log periodEnv / Prelude.fromIntegral numR
      harmonicsRD =
        fromFunction (Z :. numR :. numRFreq :. numR2Freq :. numR2Freq) $ \(Z :. r' :. rFreq' :. xFreq' :. yFreq') ->
          let rFreq = fromIntegral $ rFreq' - centerRFreq
              xFreq = fromIntegral $ xFreq' - centerR2Freq
              yFreq = fromIntegral $ yFreq' - centerR2Freq
              r = fromIntegral r' * deltaR - log periodEnv / 2
           in (1 / log periodEnv :+ 0) *
              cis
                (2 * pi / period * (x * xFreq + y * yFreq) +
                 2 * pi / log periodEnv * rFreq * r)
  harmonicsR <- computeUnboxedP harmonicsRD
  ys <-
    VS.toList . VS.map magnitude <$>
    gemmBLAS
      numR
      1
      (numRFreq * numR2Freq ^ 2)
      (VS.convert . toUnboxed $ harmonicsR)
      asteriskGaussianFreqR
  toFile def (folderPath </> "R.png") $ do
    layout_title .= printf "(%.3f , %.3f)" x y
    plot
      (line
         ""
         [ L.zip
             [ fromIntegral i * deltaR - log periodEnv / 2
             | i <- [0 .. numR - 1]
             ]
             ys
         ])
