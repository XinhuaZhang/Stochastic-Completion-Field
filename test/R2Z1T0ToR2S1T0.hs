{-# LANGUAGE FlexibleContexts #-}
module R2Z1T0ToR2S1T0 where

import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.List                 as L
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freqs = read theta0FreqsStr :: [Double]
      thetaFreqs = read thetaFreqsStr :: [Double]
      numThread = read numThreadStr :: Int
  arrR2Z1T0 <- computeR2Z1T0Array numPoint numPoint (-1) thetaFreqs theta0Freqs
    -- solveMonteCarloR2Z1T0
    --   numThread
    --   numTrail
    --   maxTrail
    --   numPoint
    --   numPoint
    --   sigma
    --   tao
    --   len
    --   theta0Freqs
    --   thetaFreqs
    --   init
  let arr = r2z1t0Tor2s1t0 numOrientation thetaFreqs arrR2Z1T0
      arr3d =
        rotate3D . R.slice arr $
        (Z :. All :. (L.length theta0Freqs - 1) :. All :. All)
      arr' =
        computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS $ arr3d
      folderPath = "output/test/R2Z1T0ToR2S1T0"
  createDirectoryIfMissing True (folderPath </> "GreensR2S1T0")
  plotImageRepaComplex (folderPath </> "GreensR2S1T0.png") . ImageRepa 8 $ arr'
  MP.mapM_
    (\i ->
       plotImageRepaComplex
         (folderPath </> "GreensR2S1T0" </> show (i + 1) L.++ ".png") .
       ImageRepa 8 .
       computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.slice arr3d $
       (Z :. All :. All :. i))
    [0 .. numOrientation - 1]
