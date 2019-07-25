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
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:numTrailStr:maxTrailStr:initOriStr:initSpeedStr:writeFlagStr:numIterationStr:thresholdStr:histFilePath:rStr:numThreadStr:_) <-
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
      r = read rStr :: Int
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2S1PointSet"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  arrG' <-
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
            r
            histFilePath
            (0, 0, 0, 0, initOri / 180 * pi, initSpeed)
  let arrG =
        computeS . R.map (\x ->  1 - x) 
        . R.traverse arrG' id $ \f idx@(Z :. k :. i :. j) ->
          let r =
                sqrt . fromIntegral $
                (i - div numPoint 2) ^ 2 + (j - div numPoint 2) ^ 2
           in if r <= 12
                 then f idx
                 else 0
  let r = 25
      numTheta = 15
      deltaTheta = (1 * pi) / numTheta
      xs = [R2S1RPPoint (-i, 0, 0, 1) | i <- [-24,-18 .. 24]]
        -- ((L.map
        --     (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
        --     [ (r * cos (k * deltaTheta) + 0, r * sin (k * deltaTheta) + 0)
        --     | k <- [0 .. numTheta]
        --     ]) -- L.++
        --  -- [R2S1RPPoint (-3, -1, 0, 1)]
        --  )
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
