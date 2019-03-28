module STCR2Z1T0ImageSparseCode where

import           Control.Monad                as M
import           Control.Monad.Parallel       as MP
import           Data.Array.Repa              as R
import           Data.Binary                  (decodeFile)
import           Data.Complex
import           Data.Conduit
import           Data.Conduit.List            as CL
import           Data.List                    as L
import           Data.Maybe
import           Data.Vector.Generic          as VG
import           Data.Vector.Storable         as VS
import           DFT.Plan
import           Filter.Gaussian
import           Filter.Pinwheel
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           Image.Transform
import           PVPFile.IO
import           STC.CompletionField
import           STC.OrientationScaleAnalysis
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initialScaleStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:histFilePath:alphaStr:pinwheelFlagStr:imagePath:reconPath:numIterationStr:writeSourceFlagStr:thresholdStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      initialScale = read initialScaleStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      alpha = read alphaStr :: Double
      pinwheelFlag = read pinwheelFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      threshold = read thresholdStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z1T0ImageSparseCode"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  imgHeader <- readPVPHeader imagePath
  reconHeader <- readPVPHeader reconPath
  let m = min (nBands imgHeader) (nBands reconHeader)
  print m
  xs <-
    runConduitRes $
    pvpFileSource imagePath .| (CL.drop (m - 1) >> CL.map id) .| CL.head
  ys <-
    runConduitRes $
    pvpFileSource reconPath .| (CL.drop (m - 1) >> CL.map id) .| CL.head
  let img' = pvpOutputData2Array . fromJust $ xs
      recon = pvpOutputData2Array . fromJust $ ys
      (Z :. cols :. rows :. _) = extent img'
  print . extent $ img'
  plotImageRepa
    (folderPath </> "recon.png")
    (ImageRepa 8 . computeS . rotate3DR $ recon)
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        solveMonteCarloR2Z1T0Radial
          numThread
          numTrail
          maxTrail
          numPoint
          numPoint
          sigma
          tao
          1
          theta0Freqs
          thetaFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length theta0Freqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z1T0 <-
    computeUnboxedP $
    computeR2Z1T0ArrayRadial
      radialArr
      numPoint
      numPoint
      1
      thetaFreqs
      theta0Freqs
  -- let arr3d =
  --       rotate3D . R.slice arrR2Z1T0 $
  --       (Z :. All :. (L.length theta0Freqs - 1) :. All :. All)
  -- MP.mapM_
  --   (\i ->
  --      plotImageRepaComplex
  --        (folderPath </> "GreensR2Z1T0_" L.++ show (i + 1) L.++ ".png") .
  --      ImageRepa 8 .
  --      computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr3d $
  --      (Z :. All :. All :. i))
  --   [0 .. (L.length thetaFreqs) - 1]
  plan <- makeR2Z1T0Plan emptyPlan arrR2Z1T0
  img'' <-
    fmap (fromUnboxed (Z :. cols :. rows)) .
    dropPixel (4 / 5) . toUnboxed . R.sumS $
    img'
  let s = R.sumAllS img'
      img = R.map (/ s) . R.sumS $ img'
      initialDistSource =
        R.traverse
          img
          (\(Z :. nx :. ny) -> (Z :. (L.length theta0Freqs) :. nx :. ny))
          (\f (Z :. _ :. i :. j) -> f (Z :. i :. j) :+ 0)
  powerMethod1'
    plan
    folderPath
    cols
    rows
    numOrientation
    thetaFreqs
    theta0Freqs
    arrR2Z1T0
    numIteration
    writeSourceFlag
    -- ""
    threshold
    initialDistSource
