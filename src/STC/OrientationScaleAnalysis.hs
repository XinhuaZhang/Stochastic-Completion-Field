{-# LANGUAGE FlexibleContexts #-}

module STC.OrientationScaleAnalysis where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           FokkerPlanck.DomainChange (r2z1Tor2s1)
import           Graphics.Gnuplot.Simple
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           STC.Convolution       
import           System.FilePath
import           Types
import           Utils.Array
import           Utils.Parallel
import Text.Printf


-- {-# INLINE analyzeOrientation #-}
-- analyzeOrientation :: Int -> [Double] -> R2T0Array -> (R2T0Array, R2T0Array)
-- analyzeOrientation numOrientation thetaFreqs arr =
--   let arrR2S1 = R.map magnitude . r2z1Tor2s1 numOrientation thetaFreqs $ arr
--       deltaTheta = 2 * pi / fromIntegral numOrientation :: Double
--       (Z :. _ :. cols :. rows) = extent arrR2S1
--       orientationArr =
--         fromListUnboxed (Z :. cols :. rows) .
--         parMap
--           rdeepseq
--           (\(i, j) ->
--              deltaTheta *
--              (fromIntegral .
--               VU.maxIndex . toUnboxed . computeS . R.slice arrR2S1 $
--               (Z :. All :. i :. j))) $
--         [(i, j) | i <- [0 .. cols - 1], j <- [0 .. rows - 1]]
--       freqArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
--       magArr =
--         R.map sqrt . sumS . R.map (\x -> (magnitude x) ^ 2) . rotate3D $ arr
--       func theta =
--         computeS .
--         R.traverse3 orientationArr freqArr magArr (\_ _ _ -> extent arr) $ \fOri fFreq fMag (Z :. k :. i :. j) ->
--           (fMag (Z :. i :. j) :+ 0) *
--           (exp $ 0 :+ (1) * fFreq (Z :. k) * (fOri (Z :. i :. j) + theta))
--    in (func 0, func 0) 


-- {-# INLINE analyzeOrientationR2Z1T0 #-}
-- analyzeOrientationR2Z1T0 ::
--      Int -> [Double] -> [Double] -> R2T0Array -> (R2Z1T0Array, R2Z1T0Array)
-- analyzeOrientationR2Z1T0 numOrientation thetaFreqs theta0Freqs arr =
--   let arrR2S1 = R.map magnitude . r2z1Tor2s1 numOrientation thetaFreqs $ arr
--       deltaTheta = 2 * pi / fromIntegral numOrientation :: Double
--       (Z :. _ :. cols :. rows) = extent arrR2S1
--       orientationArr =
--         fromListUnboxed (Z :. cols :. rows) .
--         parMap
--           rdeepseq
--           (\(i, j) ->
--              deltaTheta *
--              (fromIntegral .
--               VU.maxIndex . toUnboxed . computeS . R.slice arrR2S1 $
--               (Z :. All :. i :. j))) $
--         [(i, j) | i <- [0 .. cols - 1], j <- [0 .. rows - 1]]
--       freqArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
--       freq0Arr = fromListUnboxed (Z :. (L.length theta0Freqs)) theta0Freqs
--       magArr =
--         R.map sqrt . sumS . R.map (\x -> (magnitude x) ^ 2) . rotate3D $ arr
--       func theta =
--         computeS .
--         R.traverse4
--           orientationArr
--           freqArr
--           freq0Arr
--           magArr
--           (\_ _ _ _ ->
--              (Z :. (L.length thetaFreqs) :. (L.length theta0Freqs) :. cols :.
--               rows)) $ \fOri fFreq fFreq0 fMag (Z :. k :. l :. i :. j) ->
--           (fMag (Z :. i :. j) :+ 0) *
--           (exp $ 0 :+ (-1) * (fFreq (Z :. k) + fFreq0 (Z :. l)) * (fOri (Z :. i :. j) + theta))
--    in (func (pi / 2),func 0)


-- {-# INLINE normalizeList #-}
-- normalizeList :: (Ord e, Fractional e) => [e] -> [e]
-- normalizeList xs = L.map (/ L.maximum xs) xs

-- plotMagnitudeOrientation ::
--      (R.Source s (Complex Double))
--   => FilePath
--   -> Int
--   -> [Double]
--   -> R.Array s DIM3 (Complex Double)
--   -> (Int, Int)
--   -> IO (Int,Int)
-- plotMagnitudeOrientation folderPath numOrientationSample thetaFreqs arr (i', j') = do
--   let orientationSampleRad =
--         [ 2 * pi / fromIntegral numOrientationSample * fromIntegral i
--         | i <- [0 .. numOrientationSample - 1]
--         ]
--       orientationSampleDeg =
--         [ 360 * fromIntegral i / fromIntegral numOrientationSample
--         | i <- [0 .. numOrientationSample - 1]
--         ]
--       freqArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
--       xsFreqDomain
--         -- normalizeList $
--        =
--         parMap
--           rdeepseq
--           (\theta ->
--              magnitude .
--              R.sumAllS .
--              R.zipWith (\freq x -> x * exp (0 :+ theta * freq)) freqArr .
--              R.slice arr $
--              (Z :. All :. i :. j))
--           orientationSampleRad
--       xs =
--         normalizeList $
--         R.toList .
--         R.map magnitude .
--         r2z1Tor2s1 numOrientationSample thetaFreqs .
--         extend (Z :. All :. (1 :: Int) :. (1 :: Int)) . R.slice arr $
--         (Z :. All :. i :. j)
--       (Z :. _ :. cols :. rows) = extent arr
--       magVec =
--         VU.concat $
--         parMap
--           rdeepseq
--           (\theta ->
--              toUnboxed .
--              computeS .
--              R.map magnitude . R.sumS . rotate3D . R.traverse2 arr freqArr const $ \f1 f2 idx@(Z :. k :. _ :. _) ->
--                f1 idx * exp (0 :+ theta * f2 (Z :. k)))
--           orientationSampleRad
--       (Z :. c :. a :. b) =
--         fromIndex (Z :. (L.length thetaFreqs) :. cols :. rows) . VU.maxIndex $
--         magVec
--       maxMag = VU.maximum magVec
--       (i,j) = (a,b)
--   printf
--     "Max magnitude: %0.5f at (%d,%d) %f degree.\n"
--     maxMag
--     a
--     b
--     ((fromIntegral c :: Double) / (fromIntegral numOrientationSample) * 360)
--   plotPathsStyle
--     [ PNG (folderPath </> "Magnitude.png")
--     , Title ("Magnitude at " L.++ show (i, j))
--     ] $
--     L.zip
--         -- defaultStyle
--       --     { plotType = LinesPoints
--       --     , lineSpec = CustomStyle [LineTitle "Spatial Domain", PointType 1]
--       --     }
--       -- ,
--       [ defaultStyle
--           { plotType = Lines
--           , lineSpec =
--               CustomStyle
--                 [ LineTitle "Frequency Domain" -- , PointType 0
--                 ]
--           }
--       ]
--        -- L.zip orientationSampleDeg xs,
--       [L.zip orientationSampleDeg xsFreqDomain]
--   return (a,b)


-- plotMagnitudeOrientationSource ::
--      (R.Source s (Complex Double))
--   => DFTPlan
--   -> FilePath
--   -> Int
--   -> Int
--   -> [Double]
--   -> R2Z1T0Array
--   -> R.Array s DIM3 (Complex Double)
--   -> (Int, Int)
--   -> IO ()
-- plotMagnitudeOrientationSource plan folderPath numOrientationSample numOrientation thetaFreqs filter input (i, j) = do
--   let orientationSampleRad =
--         [ 2 * pi / fromIntegral numOrientationSample * fromIntegral i
--         | i <- [0 .. numOrientationSample - 1]
--         ]
--       orientationSampleDeg =
--         [ 360 * fromIntegral i / fromIntegral numOrientationSample
--         | i <- [0 .. numOrientationSample - 1]
--         ]
--       freqArr = fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs
--       (Z :. numThetaFreq :. cols :. rows) = extent input
--   xs <-
--     MP.mapM
--       (\theta -> do
--          initialDistF <-
--            fmap (fromUnboxed (Z :. numThetaFreq :. cols :. rows) . VU.convert) .
--            dftExecute plan (DFTPlanID DFT1DG [numThetaFreq, cols, rows] [1, 2]) .
--            VU.convert . toUnboxed . computeS . R.traverse2 input freqArr const $
--            (\f1 f2 idx@(Z :. k :. i :. j) ->
--               f1 idx * exp (0 :+ theta * (1) * f2 (Z :. k)))
--          sourceArr <- convolveR2T0 plan filter initialDistF
--          let sourceR2Z1 = R.sumS . rotateR2Z1T0Array $ sourceArr
--              mag =
--                magnitude .
--                R.sumAllS .
--                r2z1Tor2s1 numOrientation thetaFreqs .
--                extend (Z :. All :. (1 :: Int) :. (1 :: Int)) .
--                R.slice sourceR2Z1 $
--                (Z :. All :. i :. j)
--          return mag)
--       orientationSampleRad
--   plotPath
--     [ PNG (folderPath </> "SourceMagnitude.png")
--     , Title ("Source Magnitude at " L.++ show (i, j))
--     ] .
--     L.zip orientationSampleDeg . normalizeList $
--     xs
