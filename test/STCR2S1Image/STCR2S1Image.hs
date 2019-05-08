module STCR2S1Image where

import           Control.Monad           as M
import           Control.Monad.Parallel  as MP
import           Data.Array.Repa         as R
import           Data.Binary             (decodeFile)
import           Data.Complex
import           Data.List               as L
import           Data.Vector.Unboxed     as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           STC.CompletionFieldR2S1
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:numTrailStr:maxTrailStr:initOriStr:initSpeedStr:writeFlagStr:numIterationStr:thresholdStr:histFilePath:cutoffStr:imagePath:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      initOri = read initOriStr :: Double
      initSpeed = read initSpeedStr :: Double
      writeFlag = read writeFlagStr :: Bool
      numIteration = read numIterationStr :: Int
      threshold = read thresholdStr :: Double
      cutoff = read cutoffStr :: Int
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2S1Image"
  imgRepa@(ImageRepa _ img') <- readImageRepa imagePath False
  let (Z :. _ :. cols :. rows) = extent img'
      img =
        computeS . R.traverse img' id $ \f idx@(Z :. k :. i :. j) ->
          if (sqrt . fromIntegral $ (i - div cols 2) ^ 2 + (j - div rows 2) ^ 2) >
             32
            then 0
            else f idx
  createDirectoryIfMissing True folderPath
  plotImageRepa (folderPath </> "input.png") . ImageRepa 8 $ img
  flag <- doesFileExist histFilePath
  arrG <-
    if flag
      then getNormalizedHistogramArr .
           mapHistogram
             (\x ->
                let y = x :: Int
                 in fromIntegral y) <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        computeS . R.map magnitude <$>
          solveMonteCarloR2S1
            numThread
            numTrail
            numPoint
            numPoint
            numOrientation
            sigma
            tao
            cutoff
            histFilePath
            (0, 0, 0, 0, initOri / 180 * pi, initSpeed)
  powerMethod
    emptyPlan
    folderPath
    arrG
    numIteration
    writeFlag
    ""
    threshold
    (R.traverse img (const (Z :. numOrientation :. cols :. rows)) $ \f (Z :. _ :. i :. j) ->
       if f (Z :. (0 :: Int) :. i :. j) > 0
         then 1
         else 0)
    (R.traverse img (const (Z :. numOrientation :. cols :. rows)) $ \f (Z :. _ :. i :. j) ->
       if f (Z :. (0 :: Int) :. i :. j) > 0
         then 1 / fromIntegral numOrientation
         else 0)
