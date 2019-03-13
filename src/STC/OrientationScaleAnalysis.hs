{-# LANGUAGE FlexibleContexts #-}

module STC.OrientationScaleAnalysis where

import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           FokkerPlanck.DomainChange (r2z1Tor2s1)
import           Graphics.Gnuplot.Simple
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           System.FilePath
import           Types
import           Utils.Parallel
import Data.Vector.Unboxed as VU
import DFT.Plan
import Control.Monad.Parallel as MP
import STC.CompletionField
import           Utils.Array

{-# INLINE normalizeList #-}
normalizeList :: (Ord e, Fractional e) => [e] -> [e]
normalizeList xs = L.map (/ L.maximum xs) xs

plotMagnitudeOrientation ::
     (R.Source s (Complex Double))
  => FilePath
  -> Int
  -> Int
  -> [Double]
  -> R.Array s DIM3 (Complex Double)
  -> (Int, Int)
  -> IO ()
plotMagnitudeOrientation folderPath numOrientationSample numOrientation thetaFreqs arr (i, j) = do
  let orientationSampleRad =
        [ 2 * pi / fromIntegral numOrientationSample * fromIntegral i
        | i <- [0 .. numOrientationSample - 1]
        ]
      orientationSampleDeg =
        [ 360 * fromIntegral i / fromIntegral numOrientationSample
        | i <- [0 .. numOrientationSample - 1]
        ]
      freqArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
      xsFreqDomain =
        normalizeList $
        parMap
          rdeepseq
          (\theta ->
             magnitude .
             R.sumAllS .
             R.zipWith (\freq x -> x * exp (0 :+ theta * (-freq))) freqArr .
             R.slice arr $
             (Z :. All :. i :. j))
          orientationSampleRad
      xs =
        normalizeList $
        R.toList .
        R.map magnitude .
        r2z1Tor2s1 numOrientation thetaFreqs .
        extend (Z :. All :. (1 :: Int) :. (1 :: Int)) . R.slice arr $
        (Z :. All :. i :. j)
  plotPathsStyle
    [ PNG (folderPath </> "Magnitude.png")
    , Title ("Magnitude at " L.++ show (i, j))
    ] $
    L.zip
      [ defaultStyle
          { plotType = LinesPoints
          , lineSpec = CustomStyle [LineTitle "Spatial Domain", PointType 1]
          }
      , defaultStyle
          { plotType = LinesPoints
          , lineSpec = CustomStyle [LineTitle "Frequency Domain", PointType 2]
          }
      ]
      [L.zip orientationSampleDeg xs, L.zip orientationSampleDeg xsFreqDomain]


plotMagnitudeOrientationSource ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> [Double]
  -> R2Z1T0Array
  -> R.Array s DIM3 (Complex Double)
  -> (Int, Int)
  -> IO ()
plotMagnitudeOrientationSource plan folderPath numOrientationSample numOrientation thetaFreqs filter input (i, j) = do
  let orientationSampleRad =
        [ 2 * pi / fromIntegral numOrientationSample * fromIntegral i
        | i <- [0 .. numOrientationSample - 1]
        ]
      orientationSampleDeg =
        [ 360 * fromIntegral i / fromIntegral numOrientationSample
        | i <- [0 .. numOrientationSample - 1]
        ]
      freqArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
      (Z :. numThetaFreq :. cols :. rows) = extent input
  xs <-
    MP.mapM
      (\theta -> do
         initialDistF <-
           fmap (fromUnboxed (Z :. numThetaFreq :. cols :. rows) . VU.convert) .
           dftExecute plan (DFTPlanID DFT1DG [numThetaFreq, cols, rows] [1, 2]) .
           VU.convert . toUnboxed . computeS . R.traverse2 input freqArr const $
           (\f1 f2 idx@(Z :. k :. i :. j) ->
              f1 idx * exp (0 :+ theta * (1) * f2 (Z :. k)))
         sourceArr <- convolveR2T0 plan filter initialDistF
         let sourceR2Z1 = R.sumS . rotateR2Z1T0Array $ sourceArr
             mag =
               magnitude .
               R.sumAllS .
               r2z1Tor2s1 numOrientation thetaFreqs .
               extend (Z :. All :. (1 :: Int) :. (1 :: Int)) .
               R.slice sourceR2Z1 $
               (Z :. All :. i :. j)
         return mag)
      orientationSampleRad
  plotPath
    [ PNG (folderPath </> "SourceMagnitude.png")
    , Title ("Source Magnitude at " L.++ show (i, j))
    ] .
    L.zip orientationSampleDeg . normalizeList $
    xs
