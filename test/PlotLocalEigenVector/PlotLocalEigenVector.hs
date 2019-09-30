module PlotLocalEigenVector where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array
import           Utils.Parallel

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:thetaFreqsStr:scaleFreqsStr:histFileName:numThreadStr:_) <-
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
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Double
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numThread = read numThreadStr :: Int
      folderPath = "output/test/PlotLocalEigenVector"
      histFilePath = folderPath </> histFileName
      imageFolder =
        printf
          "%d_%.2f_%.2f_%.0f_%.0f_%.0f_%.0f"
          numPoint
          thetaSigma
          scaleSigma
          maxScale
          tao
          thetaFreq
          scaleFreq
  removePathForcibly (folderPath </> imageFolder)
  createDirectoryIfMissing True (folderPath </> imageFolder)
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
          thetaFreqs
          thetaFreqs
          scaleFreqs
          scaleFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scaleFreqs
             , L.length thetaFreqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  let localEigenVector =
        computeLocalEigenVector
          (ParallelParams numThread 0)
          (pinwheelHollow 8)
          (cutoff 64 radialArr)
          numPoint
          numPoint
          maxScale
          thetaFreqs
          scaleFreqs
      rotationDeg = 90
      rotatedLocalEigenVector =
        R.traverse2
          localEigenVector
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          const $ \f fFreq idx@(Z :. k :. l :. i :. j) ->
          f idx *
          ((exp (0 :+ (-1) * (fFreq (Z :. k)) * rotationDeg / 180 * pi)) --  +
           -- (exp (0 :+ (-1) * (fFreq (Z :. k)) * (rotationDeg + 90) / 180 * pi)) +
           -- (exp (0 :+ (-1) * (fFreq (Z :. k)) * (rotationDeg + 180) / 180 * pi)) +
           -- (exp (0 :+ (-1) * (fFreq (Z :. k)) * (rotationDeg + 270) / 180 * pi))
           )
  MP.mapM_
    (\(tf, sf) ->
       plotImageRepaComplex
         (folderPath </> imageFolder </> (printf "%d_%d.png" tf sf)) .
       ImageRepa 8 .
       computeS .
       R.extend (Z :. (1 :: Int) :. All :. All) . R.slice localEigenVector $
       (Z :. tf :. sf :. All :. All))
    [ (tf, sf)
    | tf <- [0 .. L.length thetaFreqs - 1]
    , sf <- [0 .. L.length scaleFreqs - 1]
    ]
  plotImageRepaComplex (folderPath </> imageFolder </> "sumFreq.png") .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.sumS . R.sumS . rotate4D . rotate4D $
    localEigenVector
  plotImageRepaComplex (folderPath </> imageFolder </> "sumFreqRotation.png") .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.sumS . R.sumS . rotate4D . rotate4D $
    rotatedLocalEigenVector
  let localEigenVectorR2S1RP =
        r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
        localEigenVector
      rotatedLocalEigenVectorR2S1RP =
        r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
        rotatedLocalEigenVector
  plotImageRepaComplex (folderPath </> imageFolder </> "sum.png") .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.sumS . R.sumS . rotate4D . rotate4D $
    localEigenVectorR2S1RP
  plotImageRepaComplex (folderPath </> imageFolder </> "sumRotation.png") .
    ImageRepa 8 .
    computeS .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.sumS . R.sumS . rotate4D . rotate4D $
    localEigenVectorR2S1RP
  MP.mapM_
    (\i ->
       plotImageRepaComplex
         (folderPath </> imageFolder </> printf "Theta%03d.png" i) .
       ImageRepa 8 .
       computeS .
       R.extend (Z :. (1 :: Int) :. All :. All) . R.slice localEigenVectorR2S1RP $
       (Z :. i :. (0 :: Int) :. All :. All))
    [0 .. numOrientation - 1]
