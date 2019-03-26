module EquivarianceDiagram where

import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           DFT.Plan
import           Filter.Gaussian
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath

main = do
  (inputPath:sigmaStr:_) <- getArgs
  (ImageRepa _ img) <- readImageRepa inputPath True
  let (Z :. channels :. cols :. rows) = extent img
      sigma = read sigmaStr :: Double
      folderPath = "output/test/EquivarianceDiagram"
  createDirectoryIfMissing True folderPath
  lock <- getFFTWLock
  (plan1, imgF) <-
    dft1dGPlan lock emptyPlan [channels, cols, rows] [1, 2] .
    VU.convert . VU.map (:+ 0) . toUnboxed $
    img
  (plan2, _) <- idft1dGPlan lock plan1 [channels, cols, rows] [1, 2] imgF
  (plan, gaussianFilterF) <-
    gaussian2DFilter plan2 (Gaussian2DParams sigma rows cols)
  convolvedImg <- convolveGaussian2D plan gaussianFilterF . R.map (:+ 0) $ img
  plotImageRepa (folderPath </> (takeBaseName inputPath) L.++ "_out.png") .
    ImageRepa 8 . computeS . R.map magnitude $
    convolvedImg
