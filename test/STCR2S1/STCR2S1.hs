module STCR2S1 where

import           Control.Monad             as M
import           Control.Monad.Parallel    as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           STC.CompletionFieldR2S1
import           System.Directory
import           System.Environment
import           System.FilePath
import           Types
import           Utils.Array
import STC.Utils

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:numTrailStr:maxTrailStr:initDistStr:initOriStr:initSpeedStr:rStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      initDist = read initDistStr :: [R2S1RPPoint]
      initOri = read initOriStr :: Double
      initSpeed = read initSpeedStr :: Double
      r = read rStr :: Int
      numThread = read numThreadStr :: Int
      sourceDist =
        computeInitialDistribution numPoint numPoint numOrientation . L.take 1 $
        initDist
      sinkDist =
        computeInitialDistribution numPoint numPoint numOrientation . L.drop 1 $
        initDist
      dist =
        computeInitialDistribution numPoint numPoint numOrientation initDist
      folderPath = "output/test/STCR2S1"
  createDirectoryIfMissing True folderPath
  -- Compute the Green's function
  arrG <-
    solveMonteCarloR2S1
      numThread
      numTrail
      numPoint
      numPoint
      numOrientation
      sigma
      tao
      r
      ""
  plan <- makeR2S1Plan emptyPlan arrG
  -- Source Field
  source <- shareWeightST plan sourceDist arrG
  plotImageRepa
    (folderPath </> "Source.png")
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D . reduceContrast 50$
     source)
  -- Sink Field
  sink <- shareWeightST plan sinkDist . rotateST arrG $ div numOrientation 2
  plotImageRepa
    (folderPath </> "Sink.png")
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D . reduceContrast 50$
     sink)
  -- Completion Field
  let completion = R.zipWith (*) source . timeReversal $ sink
  plotImageRepa
    (folderPath </> "Completion.png")
    (ImageRepa 8 .
     computeS . R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D . reduceContrast 50 $
     completion)
