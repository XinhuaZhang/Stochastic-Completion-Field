module OrientationAnalysisShape where

import           Control.Monad           as M
import           Data.Array.Repa         as R
import           Data.Complex
import           Data.List               as L
import           DFT.Plan
import           Filter.Pinwheel
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC
import           STC.OrientationScaleAnalysis
import           System.Directory
import           System.Environment
import           System.FilePath

main = do
  args@(numPointStr:theta0FreqsStr:alphaStr:idxStr:numOrientationSampleStr:shapeStr:whiteBackgroundFlagStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      alpha = read alphaStr :: Double
      idx = read idxStr :: (Int, Int)
      numOrientationSample = read numOrientationSampleStr :: Int
      whiteBackgroundFlag = read whiteBackgroundFlagStr :: Bool
      numThread = read numThreadStr :: Int
      pinwheelParams =
        PinwheelParams numPoint numPoint alpha (exp 1) theta0Freqs [0]
      shape2D = read shapeStr :: Shape2D
      shapeArr' =
        makeShape2DContrast numPoint numPoint (ContrastArea whiteBackgroundFlag shape2D)
      folderPath = "output/test/OrientationAnalysisShape"
  createDirectoryIfMissing True folderPath
  shapeArr <- computeUnboxedP shapeArr'
  
  (plan1, imgF) <-
    makeImagePlan emptyPlan .
    computeS .
    R.map fromIntegral .
    R.backpermute
      (Z :. (3 :: Int) :. numPoint :. numPoint)
      (\(Z :. _ :. i :. j) -> (Z :. i :. j)) $
    shapeArr
  (plan, filterF, filterPIF) <- pinwheelFilter plan1 pinwheelParams
  convolvedImg <-
    fmap (\x -> R.slice x (Z :. All :. (0 :: Int) :. All :. All)) $
    convolvePinwheel plan filterF (R.slice imgF (Z :. (0 :: Int) :. All :. All))
  idx' <- plotMagnitudeOrientation
            folderPath
            numOrientationSample
            theta0Freqs
            convolvedImg
            idx
  plotImageRepa (folderPath </> "Shape.png") .
    ImageRepa 8 .
    computeS .
    R.traverse
      shapeArr
      (\(Z :. cols :. rows) -> (Z :. (3 :: Int) :. cols :. rows)) $
    (\f (Z :. k :. i :. j) ->
       if (i, j) == idx'
         then if k == 0
                then 255
                else 0
         else fromIntegral $ f (Z :. i :. j))
