{-# LANGUAGE FlexibleContexts #-}
module PlotGreensFunctionSTS where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Unboxed       as VU
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array

main = do
  args@(numPointStr:spatialFreqStr:numOrientationStr:thetaSigmaStr:numScaleStr:scaleSigmaStr:maxScaleStr:taoStr:initStr:numTrailStr:maxTrailStr:thetaFreqsStr:scaleFreqsStr:histFileName:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      spatialFreq = read spatialFreqStr :: Double
      spatialFreqs = [-spatialFreq .. spatialFreq]
      numOrientation = read numOrientationStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      numScale = read numScaleStr :: Int
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      init@(initX, initY, initTheta, initScale) =
        read initStr :: (Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      thetaFreq = read thetaFreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numThread = read numThreadStr :: Int
      folderPath = "output/test/PlotGreensFunctionSTS"
      histPath = folderPath </> histFileName
      initArr =
        traverse4
          (fromListUnboxed (Z :. (L.length spatialFreqs)) spatialFreqs)
          (fromListUnboxed (Z :. (L.length spatialFreqs)) spatialFreqs)
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          (fromListUnboxed (Z :. (L.length scaleFreqs)) scaleFreqs)
          (\(Z :. a) (Z :. b) (Z :. c) (Z :. d) ->
             (Z :. a :. b :. c :. d :. c :. d)) $ \fx fy ft fs (Z :. x :. y :. _ :. _ :. t :. s) ->
          exp
            (0 :+
             (-1) *
             (ft (Z :. t) * (initTheta / 180 * pi) +
              2 * pi *
              ((fx (Z :. x) * initX + fy (Z :. y) * initY) /
               (fromIntegral numPoint) +
               fs (Z :. s) * initScale / log maxScale)))
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histPath
  arr <-
    if flag
      then fmap
             (fromListUnboxed
                (Z :. (L.length spatialFreqs) :. (L.length spatialFreqs) :.
                 (L.length thetaFreqs) :.
                 (L.length scaleFreqs) :.
                 (L.length thetaFreqs) :.
                 (L.length scaleFreqs))) .
           decodeFile $
           histPath
      else runMonteCarloSTS
             numThread
             numTrail
             maxTrail
             numPoint
             numPoint
             thetaSigma
             scaleSigma
             maxScale
             tao
             spatialFreqs
             spatialFreqs
             thetaFreqs
             scaleFreqs
             histPath $
           VU.replicate
             (L.length scaleFreqs * L.length thetaFreqs * L.length scaleFreqs *
              L.length thetaFreqs *
              L.length spatialFreqs *
              L.length spatialFreqs)
             0
  transformedArr <- sumP . sumS $ R.zipWith (*) arr initArr
  -- arrR2S1RP <-
  --   stsTor2s1rp
  --     numPoint
  --     (fromIntegral numPoint)
  --     spatialFreqs
  --     numPoint
  --     (fromIntegral numPoint)
  --     spatialFreqs
  --     numOrientation
  --     thetaFreqs
  --     numScale
  --     scaleFreqs
  --     maxScale
  --     transformedArr
  -- plotImageRepaComplex (folderPath </> "Greens.png") .
  --   ImageRepa 8 .
  --   computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . sumS $
  --   arrR2S1RP
  arrMag <-
    stsTor2'
      numPoint
      (fromIntegral numPoint)
      spatialFreqs
      numPoint
      (fromIntegral numPoint)
      spatialFreqs
      transformedArr
  plotImageRepa (folderPath </> "GreensMag.png") .
    ImageRepa 8 .
    computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.map sqrt $
    arrMag
