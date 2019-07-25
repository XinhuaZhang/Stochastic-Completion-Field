module RadialPhase where

import           Control.Monad           as M
import           Control.Monad.Parallel  as MP
import           Data.Array.Repa         as R
import           Data.Binary
import           Data.Complex
import           Data.List               as L
import           Data.Vector.Unboxed     as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Graphics.Gnuplot.Simple
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:maxScaleStr1:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:alphaStr:deltaThetaStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      maxScale1 = read maxScaleStr1 :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      deltaTheta = read deltaThetaStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/RadialPhase"
  createDirectoryIfMissing True folderPath
  createDirectoryIfMissing True (folderPath </> "Pinwheel")
  createDirectoryIfMissing True (folderPath </> "Magnitude")
  radialArr <-
    solveMonteCarloR2Z2T0S0Radial'
      numThread
      numTrail
      maxTrail
      numPoint
      numPoint
      thetaSigma
      scaleSigma
      maxScale
      tao
      theta0Freqs
      thetaFreqs
      scale0Freqs
      scaleFreqs
      histFilePath
      (emptyHistogram
         [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
         , L.length scale0Freqs
         , L.length theta0Freqs
         , L.length scaleFreqs
         , L.length thetaFreqs
         ]
         0)
  let arr3d =
        R.slice
          radialArr
          (Z :. (L.length thetaFreqs - 1) :. All :. (L.length theta0Freqs - 1) :.
           All :.
           All)
      (Z :. _ :. _ :. _ :. _ :. len) = extent radialArr
      arrR2Z2T0S0 =
        computeR2Z2T0S0ArrayRadial'
          numPoint
          numPoint
          1
          maxScale
          thetaFreqs
          scaleFreqs
          theta0Freqs
          scale0Freqs
      arr3d' =
        R.slice
          arrR2Z2T0S0
          (Z :. (L.length thetaFreqs - 1) :. All :. (L.length theta0Freqs - 1) :.
           All :.
           (0 :: Int) :.
           All)
      (Z :. _ :. _ :. _ :. _ :. _ :. len') = extent arrR2Z2T0S0
  print len
  print len'
  MP.mapM_
    (\((i, sf), (j, s0f)) ->
       plotPathStyle
         [ PNG (folderPath </> printf "%d_%d.png" i j)
         , Title (printf "sf: %.0f, s0f: %.0f" sf s0f)
         -- , XRange (0, log (fromIntegral len))
         , YRange (0, 2 * pi)
         ]
         (defaultStyle
            {plotType = LinesPoints, lineSpec = CustomStyle [LineTitle ""]}) .
       L.zip (L.map log [1 .. fromIntegral len - 1]) .
       L.tail . R.toList . R.map (\x -> (phase x) + pi) . R.slice arr3d $
       (Z :. i :. j :. All))
    [ (a, b)
    | a <- L.zip [0 .. (L.length scaleFreqs) - 1] scaleFreqs
    , b <- L.zip [0 .. (L.length scale0Freqs) - 1] scale0Freqs
    ]
  let maxMag = VU.maximum . toUnboxed . computeS . R.map magnitude $ arr3d
  MP.mapM_
    (\((i, sf), (j, s0f)) ->
       plotPathStyle
         [ PNG (folderPath </> "Magnitude" </> printf "%d_%d.png" i j)
         , Title (printf "sf: %.0f, s0f: %.0f" sf s0f)
         -- , XRange (0, log (fromIntegral len))
         -- , YRange (0, maxMag)
         ]
         (defaultStyle
            {plotType = LinesPoints, lineSpec = CustomStyle [LineTitle ""]}) .
       L.zip ( [0 .. fromIntegral len - 1]) .
       -- L.tail .
      R.toList . R.map magnitude . R.slice arr3d $
       (Z :. i :. j :. All))
    [ (a, b)
    | a <- L.zip [0 .. (L.length scaleFreqs) - 1] scaleFreqs
    , b <- L.zip [0 .. (L.length scale0Freqs) - 1] scale0Freqs
    ]
  MP.mapM_
    (\((i, sf), (j, s0f)) ->
       plotPathStyle
         [ PNG (folderPath </> "Pinwheel" </> printf "%d_%d.png" i j)
         , Title (printf "sf: %.0f, s0f: %.0f" sf s0f)
         -- , XRange (0, log (fromIntegral len))
         , YRange (0, 2 * pi)
         ]
         (defaultStyle
            {plotType = LinesPoints, lineSpec = CustomStyle [LineTitle ""]}) .
       L.zip (L.map log [1 .. fromIntegral len - 1]) .
       L.map
         (\r ->
            pi +
            (phase $
             exp
               (0 :+
                (s0f *
                 (log r * 2 * pi / log maxScale1 + (deltaTheta / 180 * pi)) +
                 sf * (log r * 2 * pi / log maxScale1 + (deltaTheta / 180 * pi)))))) $
       [1 .. fromIntegral len - 1])
    [ (a, b)
    | a <- L.zip [0 .. (L.length scaleFreqs) - 1] scaleFreqs
    , b <- L.zip [0 .. (L.length scale0Freqs) - 1] scale0Freqs
    ]
