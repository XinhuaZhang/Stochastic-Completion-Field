module STCR2Z2T0S0EdgeBinary where

import           Control.Arrow
import           Control.Monad
import           Data.Array.Repa         as R
import           Data.Binary             (encodeFile, decodeFile)
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.Edge
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
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:thetaFreqsStr:scaleFreqsStr:hollowRadiusStr:cutoffRadiusStr:histFilePath:filterFileFolder:numIterationStr:writeSourceFlagStr:saveEdgeDataFlagStr:loadEdgeDataFlagStr:edgeFilePath:numNoisePointStr:scaleFactorStr:useFFTWWisdomFlagStr:fftwWisdomFileName:batchSizeStr:numThreadStr:_) <-
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
      hollowRadius = read hollowRadiusStr :: Double
      cutoffRadius = read cutoffRadiusStr :: Double
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      saveEdgeDataFlag = read saveEdgeDataFlagStr :: Bool
      loadEdgeDataFlag = read loadEdgeDataFlagStr :: Bool
      numNoisePoint = read numNoisePointStr :: Int
      scaleFactor = read scaleFactorStr :: Double
      useFFTWWisdomFlag = read useFFTWWisdomFlagStr :: Bool
      batchSize = read batchSizeStr :: Int
      numThread = read numThreadStr :: Int
      parallelParams = ParallelParams numThread batchSize
      folderPath = "output/test/STCR2Z2T0S0EdgeBinary"
      fftwWisdomFilePath = folderPath </> fftwWisdomFileName
      filterFileName =
        printf
          "Filter_%.0f_%.0f_%s"
          hollowRadius
          cutoffRadius
          (takeFileName histFilePath)
      filterFilePath = filterFileFolder </> filterFileName
      biasFilePath =
        printf
          "%s/Bias_%s_%d.dat"
          filterFileFolder
          (takeBaseName edgeFilePath)
          numNoisePoint
      eigenVecFilePath =
        printf
          "%s/EigenVec_%s_%d.dat"
          filterFileFolder
          (takeBaseName edgeFilePath)
          numNoisePoint
  createDirectoryIfMissing True folderPath
  createDirectoryIfMissing True filterFileFolder
  plan <-
    makePlanBinary
      emptyPlan
      useFFTWWisdomFlag
      fftwWisdomFilePath
      (L.length thetaFreqs)
      (L.length scaleFreqs)
      numPoint
      numPoint
  filterFlag <- doesFileExist filterFilePath
  flag <-
    if filterFlag
      then do
        size <- getFileSize filterFilePath
        return $
          if size == 0
            then False
            else True
      else return False
  unless
    flag
    (do histFlag <- doesFileExist histFilePath
        radialArr <-
          if histFlag
            then R.map magnitude . getNormalizedHistogramArr <$>
                 decodeFile histFilePath
            else do
              putStrLn
                "Couldn't find a Green's function data. Start simulation..."
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
        putStrLn "Write filter to disk..."
        writeDFTPinwheel
          parallelParams
          plan
          radialArr
          hollowRadius
          cutoffRadius
          thetaFreqs
          scaleFreqs
          maxScale
          (numPoint, numPoint)
          filterFilePath)
  (bias, eigenVec) <-
    if loadEdgeDataFlag
      then do
        bias <- readRepaArray biasFilePath
        eigenVec <- readRepaArray eigenVecFilePath
        return (bias, eigenVec)
      else do
        xs <- parseEdgeFile edgeFilePath
        randomPonintSet <-
          generateRandomPointSet numNoisePoint numPoint numPoint
        let (centerX, centerY) =
              join (***) (\x -> div x . L.length $ xs) .
              L.foldl'
                (\(a, b) (R2S1RPPoint (c, d, _, _)) -> (a + c, b + d))
                (0, 0) $
              xs
            ys =
              L.map
                (\(R2S1RPPoint (a, b, c, d)) ->
                   (R2S1RPPoint
                      ( round $ (fromIntegral $ a - centerX) / scaleFactor
                      , round $ (fromIntegral $ b - centerX) / scaleFactor
                      , c
                      , d)))
                xs
            zs =
              L.map
                (\(R2S1RPPoint (a, b, c, d)) ->
                   (R2S1RPPoint (a - center numPoint, b - center numPoint, c, d)))
                randomPonintSet
            points = ys L.++ zs
        let bias =
              computeBiasR2T0S0 numPoint numPoint thetaFreqs scaleFreqs points
            eigenVec =
              computeInitialEigenVectorBinary
                numPoint
                numPoint
                thetaFreqs
                scaleFreqs
                points
        when saveEdgeDataFlag (writeRepaArray biasFilePath bias)
        return (bias, eigenVec)
  powerMethodBinary
    parallelParams
    plan
    folderPath
    numPoint
    numPoint
    numOrientation
    thetaFreqs
    numScale
    scaleFreqs
    maxScale
    filterFilePath
    numIteration
    writeSourceFlag
    (printf
       "_%.2f_%.2f_%d_%d_%d_%d"
       thetaSigma
       scaleSigma
       (round maxScale :: Int)
       (round tao :: Int)
       (round thetaFreq :: Int)
       (round scaleFreq :: Int))
    saveEdgeDataFlag
    eigenVecFilePath
    bias
    eigenVec
