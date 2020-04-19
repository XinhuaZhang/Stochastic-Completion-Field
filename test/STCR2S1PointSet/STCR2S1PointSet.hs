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
import STC.Shape
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array
import FokkerPlanck.Histogram
import FokkerPlanck.GreensFunction

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:numTrailStr:maxTrailStr:initOriStr:initSpeedStr:writeFlagStr:numIterationStr:thresholdStr:histFilePath:rStr:deltaStr:stdStr:numThreadStr:_) <-
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
      r = read rStr :: Double
      delta = read deltaStr :: Double
      std = read stdStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2S1PointSet"
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  -- arrG' <-
  --   if flag
  --     then getNormalizedHistogramArr <$> decodeFile histFilePath
  --     else do
  --       putStrLn "Couldn't find a Green's function data. Start simulation..."
  --       solveMonteCarloR2S1
  --         numThread
  --         numTrail
  --         maxTrail
  --         numPoint
  --         numPoint
  --         numOrientation
  --         sigma
  --         tao
  --         ((sqrt 2) * fromIntegral numPoint)
  --         initSpeed
  --         histFilePath
  --         -- (0, 0, 0, 0, initOri / 180 * pi, initSpeed)
  let arrG' =
        sampleR2S1 numPoint numPoint numOrientation delta initSpeed sigma tao
      arrG =
        computeS . R.traverse arrG' id $ \f idx@(Z :. k :. i :. j) ->
          let r' =
                sqrt . fromIntegral $
                (i - div numPoint 2) ^ 2 + (j - div numPoint 2) ^ 2
          in if r' > 91 -- || r' <= 2
              then 0
              else f idx
      -- r = 16
  let numTheta = 8
      deltaTheta = (2 * pi) / numTheta
      xs -- [R2S1RPPoint (-i, 0, 0, 1) | i <- [-24,-18 .. 24]]
       =
        ((L.map
            (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
            [ (r * cos (k * deltaTheta) + 0, r * sin (k * deltaTheta) + 0)
            | k <- [0 .. numTheta - 1]
            ]) -- L.++
         -- [R2S1RPPoint (-3, -1, 0, 1)]
         )
      bias = computeBias numPoint numPoint numOrientation xs
      -- initialEigenVec =
      --   computeInitialEigenVec numPoint numPoint numOrientation xs
      circle = Points (0, 0) 2 Circle {circleNum = 8, circleRadius = round r}
      circlePoints = makeShape2D $ circle
      ys = getShape2DIndexListGaussian 5 2 $ circlePoints
      initialEigenVec =
        computeInitialEigenVecGaussian numPoint numPoint numOrientation ys
  print xs
  print circlePoints
  -- powerMethod
  --   emptyPlan
  --   folderPath
  --   arrG
  --   numIteration
  --   writeFlag
  --   ""
  --   threshold
  --   bias
  --   initialEigenVec
  computeContourR2S1Analytic
    folderPath
    arrG
    numOrientation
    numPoint
    numPoint
    sigma
    tao
    initSpeed
    std
    delta
    bias
    circlePoints
  -- computeContourR2S1Tangent
  --   folderPath
  --   arrG
  --   numOrientation
  --   numPoint
  --   numPoint
  --   delta
  --   circlePoints
