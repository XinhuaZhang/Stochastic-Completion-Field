module STCR2S1PointSet where

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
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:numTrailStr:maxTrailStr:initOriStr:initSpeedStr:writeFlagStr:numIterationStr:thresholdStr:histFilePath:numThreadStr:_) <-
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
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2S1PointSet"
  createDirectoryIfMissing True folderPath
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
            histFilePath
            (0, 0, 0, 0, initOri / 180 * pi, initSpeed)
  let r = 20
      numTheta = 8
      deltaTheta = (2 * pi) / numTheta
      xs =
        ((L.map
            (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
            [ (r * cos (k * deltaTheta) + 6, r * sin (k * deltaTheta) + 2)
            | k <- [0 .. numTheta - 1]
            ]) -- L.++
         -- [R2S1RPPoint (-3, -1, 0, 1)]
        )
      bias = computeBias numPoint numPoint numOrientation xs
      initialEigenVec =
        computeInitialEigenVec numPoint numPoint numOrientation xs
  powerMethod
    emptyPlan
    folderPath
    arrG
    numIteration
    writeFlag
    ""
    threshold
    bias
    initialEigenVec
