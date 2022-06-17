{-# LANGUAGE FlexibleContexts #-}
module PlotGreensFunctionR2S1RP where
import           Control.Monad           as M
import           Control.Monad.Parallel  as MP
import           Data.Array.Repa         as R
import           Data.List               as L
import           FokkerPlanck.MonteCarlo
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
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
  (numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:lenStr:initStr:numTrailStr:numThreadStr:_) <-
    getArgs
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      len = read lenStr :: Int
      init = read initStr :: (Double, Double, Double, Double, Double, Double)
      numTrail = read numTrailStr :: Int
      numThread = read numThreadStr :: Int
  arr'' <-
    solveMonteCarloR2S1RP
      numThread
      numTrail
      numPoint
      numPoint
      numOrientation
      numScale
      thetaSigma
      scaleSigma
      maxScale
      tao
      len
      init
  let arr = rotate4D . rotate4D $  arr''
      arr' =
        computeS .
        reduceContrast 10 .
        R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . R.sumS $
        arr
      folderPath = "output/app/PlotGreensFunctionR2S1RP"
  createDirectoryIfMissing True (folderPath </> "GreensR2S1RP")
  plotImageRepa (folderPath </> "GreensR2S1RP.png") . ImageRepa 8 $ arr'
  MP.mapM_
    (\i ->
       plotImageRepa
         (folderPath </> "GreensR2S1RP/" L.++ (show $ i + 1) L.++ ".png") .
       ImageRepa 8 .
       computeS .
       reduceContrast 50 .
       R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . R.slice arr $
       (Z :. All :. All :. All :. i))
    [0 .. numScale - 1]
