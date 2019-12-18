{-# LANGUAGE FlexibleContexts #-}
module PlotGreensFunctionR2S1 where

import           Data.Array.Repa         as R
import           Data.Complex
import           Data.List               as L
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array

{-# INLINE reduceContrast #-}
reduceContrast
  :: (R.Source s Double)
  => Int -> Array s DIM3 Double -> Array D DIM3 Double
reduceContrast idx arr =
  let x = L.head . L.drop idx . L.reverse . L.sort . R.toList $ arr
  in R.map
       (\y ->
          if y >= x
            then x
            else y)
       arr

main = do
  args@(numPointStr:numOrientationStr:sigmaStr:taoStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      sigma = read sigmaStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init@ (_,_,_,_,_,s0) = read initStr :: (Double, Double, Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      numThread = read numThreadStr :: Int
  arr <-
    solveMonteCarloR2S1
      numThread
      numTrail
      numPoint
      numPoint
      numOrientation
      sigma
      tao
      numPoint
      ""
      init
  let arr' =
        computeS .
        -- reduceContrast 10 .
        R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . R.map magnitude . rotate3D $
        arr
      folderPath = "output/app/PlotGreensFunctionR2S1"
  createDirectoryIfMissing True folderPath
  plotImageRepa (folderPath </> printf "GreensR2S1_%.0f.png" s0) . ImageRepa 8 $ arr'
