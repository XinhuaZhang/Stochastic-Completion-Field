module Deconvolution where

import           Control.Monad             as M
import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Graphics.Gnuplot.Simple
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:initDistStr:histFilePath:alphaStr:pinwheelFlagStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      initDist = read initDistStr :: [R2S1RPPoint]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/Deconvolution"
      (R2S1RPPoint (_, _, _, s)) = L.head initDist
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        solveMonteCarloR2Z2T0S0Radial
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
  radialArrSink <-
    R.map magnitude . getNormalizedHistogramArr <$> decodeFile histFilePath
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      radialArr
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  arrR2Z2T0S0Sink <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (cutoff 24 radialArrSink)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  sourceDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      theta0Freqs
      scale0Freqs
      maxScale
      initDist
  sourceDistArrR2S1RP <-
    r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $ sourceDistArr
  let xIndex =
        [ (fromIntegral i) * 360 / (fromIntegral numOrientation)
        | i <- [0 .. numOrientation - 1]
        ] :: [Double]
      inputStyle =
        L.zipWith
          (\(R2S1RPPoint (x', y', _, _)) i ->
             let x = x' + center numPoint
                 y = y' + center numPoint
              in ( defaultStyle
                     { plotType = LinesPoints
                     , lineSpec =
                         CustomStyle
                           [LineTitle (printf "(%d,%d)" x y), PointType i]
                     }
                 , L.zip xIndex .
                   R.toList . R.slice (R.map magnitude sourceDistArrR2S1RP) $
                   (Z :. All :. (0 :: Int) :. x :. y)))
          initDist
          [1 ..]
  plotPathsStyle [PNG (folderPath </> "Input.png"), Title "Input"] inputStyle
  arrR2Z2T0S0F <-
    dftR2Z2T0S0 plan .
    computeS .
    makeFilter2D .
    R.traverse arrR2Z2T0S0
      -- (R.traverse (arrR2Z2T0S0) id $ \f idx@(Z :. _ :. _ :. _ :. _ :. a :. b) ->
      --    if a == center numPoint && b == center numPoint
      --      then 0
      --      else f idx)
      id $ \f (Z :. tf' :. sf' :. t0f :. s0f :. i :. j) ->
      let idx = (Z :. tf' :. sf' :. t0f :. s0f :. i :. j)
       in if i == center numPoint && j == center numPoint &&
             tf' == t0f
                    -- tf' == div (L.length thetaFreqs) 2 &&
                    -- sf' == div (L.length scaleFreqs) 2 &&
                    -- t0f == div (L.length theta0Freqs) 2 &&
                    -- s0f == div (L.length scale0Freqs) 2 
            then f idx + (0.0 :+ 0)
            else f idx
  -- let arr' =
  --       r2z2t0s0Tor2s1rps1rp
  --         numOrientation
  --         thetaFreqs
  --         theta0Freqs
  --         numScale
  --         scaleFreqs
  --         scale0Freqs
  --         maxScale $
  --       arrR2Z2T0S0Sink
  --     arr'' =
  --       r2s1rps1rpTor2z2t0s0
  --         numOrientation
  --         thetaFreqs
  --         theta0Freqs
  --         numScale
  --         scaleFreqs
  --         scale0Freqs
  --         maxScale .
  --       computeUnboxedS .
  --       R.map
  --         (\x ->
  --            let (m, p) = polar x
  --             in if m == 0
  --                  then 0
  --                  else mkPolar (1 / m) p) $
  --       arr'
  arrR2Z2T0S0SinkF <-
    dftR2Z2T0S0 plan .
    computeS .
    makeFilter2D . computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs
    $
    arrR2Z2T0S0Sink
  --   R.traverse
  --     (-- computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs $
  --      arrR2Z2T0S0)
  --     id $ \f idx@(Z :. a :. b :. c :. d :. i :. j) ->
  --     let (m, p) =
  --           polar $
  --           conjugate
  --             (f (Z :. (L.length thetaFreqs - 1 - a) :.
  --                 (L.length scaleFreqs - 1 - b) :.
  --                 (L.length theta0Freqs - 1 - c) :.
  --                 (L.length scale0Freqs - 1 - d) :.
  --                 i :.
  --                 j))
  --      in if a == div (L.length thetaFreqs) 2 &&
  --            b == div (L.length scaleFreqs) 2 &&
  --            c == div (L.length theta0Freqs) 2 &&
  --            d == div (L.length scale0Freqs) 2
  --           then 1 + conjugate
  --                       (f (Z :. (L.length thetaFreqs - 1 - a) :.
  --                           (L.length scaleFreqs - 1 - b) :.
  --                           (L.length theta0Freqs - 1 - c) :.
  --                           (L.length scale0Freqs - 1 - d) :.
  --                           i :.
  --                           j))
  --           else conjugate
  --                  (f (Z :. (L.length thetaFreqs - 1 - a) :.
  --                      (L.length scaleFreqs - 1 - b) :.
  --                      (L.length theta0Freqs - 1 - c) :.
  --                      (L.length scale0Freqs - 1 - d) :.
  --                      i :.
  --                      j))
  -- Source field
  sourceArr <- convolveR2T0S0 plan arrR2Z2T0S0F sourceDistArr
  sourceR2Z2' <- R.sumP . R.sumS . rotateR2Z2T0S0Array $ sourceArr
  sourceR2S1RP <-
    r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $ sourceR2Z2'
  let sourceR2Z2 -- = sourceR2Z2'
       =
        R.traverse2 sourceR2Z2' sourceDistArr const $ \f1 f2 idx@(Z :. a :. b :. c :. d) ->
          let x = div (L.length theta0Freqs) 2
              y = div (L.length scale0Freqs) 2
           in if f2 (Z :. x :. y :. c :. d) == 0
                then 0
                else f1 idx
      sourceStyle =
        L.zipWith
          (\(R2S1RPPoint (x', y', _, _)) i ->
             let x = x' + center numPoint
                 y = y' + center numPoint
              in ( defaultStyle
                     { plotType = LinesPoints
                     , lineSpec =
                         CustomStyle
                           [LineTitle (printf "(%d,%d)" x y), PointType i]
                     }
                 , L.zip xIndex .
                   R.toList . R.slice (R.map magnitude sourceR2S1RP) $
                   (Z :. All :. (0 :: Int) :. x :. y)))
          initDist
          [1 ..]
  plotPathsStyle [PNG (folderPath </> "Conv.png"), Title "Conv"] sourceStyle
  sourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    sourceR2Z2
  plotImageRepaComplex (folderPath </> "Source.png") . ImageRepa 8 $ sourceField
  -- Deconvolution
  deconvSourceArr <- convolveR2T0S0 plan arrR2Z2T0S0SinkF . computeUnboxedS $ sourceR2Z2
  deconvSourceR2Z2 <-
    R.sumP . R.sumS . rotateR2Z2T0S0Array $ deconvSourceArr
  deconvSourceField <-
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    deconvSourceR2Z2
  plotImageRepaComplex (folderPath </> "DeconvSource.png") . ImageRepa 8 $
    deconvSourceField
  deconvSourceR2S1RP <-
    r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
    deconvSourceR2Z2
  let outputStyle =
        L.zipWith
          (\(R2S1RPPoint (x', y', _, _)) i ->
             let x = x' + center numPoint
                 y = y' + center numPoint
              in ( defaultStyle
                     { plotType = LinesPoints
                     , lineSpec =
                         CustomStyle
                           [LineTitle (printf "(%d,%d)" x y), PointType i]
                     }
                 , L.zip xIndex .
                   R.toList . R.slice (R.map magnitude deconvSourceR2S1RP) $
                   (Z :. All :. (0 :: Int) :. x :. y)
                     -- (R.sumS .
                   --    R.backpermute
                   --      (Z :. (L.length theta0Freqs) :. numPoint :. numPoint :.
                   --       (L.length scale0Freqs))
                   --      (\(Z :. a :. b :. c :. d) -> (Z :. a :. d :. b :. c)) .
                   --    R.map magnitude $
                   --    deconvSourceR2S1RP) $
                   -- (Z :. All :. x :. y)
                  ))
          initDist
          [1 ..]
  plotPathsStyle [PNG (folderPath </> "Output.png"), Title "Output"] outputStyle
