{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
module Filter.Gaussian where

import           Control.Monad       (when)
import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List           as L (length, reverse, (++))
import           Data.Vector.Generic as VG (convert)
import           DFT.Plan
import           Filter.Utils

data Gaussian2DParams = Gaussian2DParams
  { getGaussianFilterSigma :: !Double
  , getGaussianFilterRows  :: !Int
  , getGaussianFilterCols  :: !Int
  } deriving (Show)

{-# INLINE gaussian2D #-}
gaussian2D
  :: (Floating a)
  => a -> Int -> Int -> a
gaussian2D sd i j =
  1 / ((2 * pi) * sd * sd) * exp (-r / (2 * (sd ^ (2 :: Int))))
  where
    r = fromIntegral (i * i + j * j)

gaussian2DFilter ::
     DFTPlan -> Gaussian2DParams -> IO (DFTPlan, Array U DIM2 (Complex Double))
gaussian2DFilter plan (Gaussian2DParams sigma rows cols) = do
  let filterVec =
        VG.convert .
        toUnboxed .
        computeS .
        R.map (:+ 0) . makeFilter . R.fromFunction (Z :. cols :. rows) $ \(Z :. i :. j) ->
          gaussian2D sigma (i - (div cols 2)) (j - (div rows 2))
  lock <- getFFTWLock
  (plan1, filterVecF) <- dft1dGPlan lock plan [cols, rows] [0, 1] filterVec
  (newPlan, _) <- idft1dGPlan lock plan1 [cols, rows] [0, 1] filterVecF
  return (newPlan, fromUnboxed (Z :. cols :. rows) . VG.convert $ filterVecF)

{-# INLINE convolveGaussian2D #-}
convolveGaussian2D ::
     (Source r (Complex Double), Shape sh)
  => DFTPlan
  -> Array U DIM2 (Complex Double)
  -> Array r (sh :. Int :. Int) (Complex Double)
  -> IO (Array U (sh :. Int :. Int) (Complex Double))
convolveGaussian2D plan filterF input = do
  let len = rank . extent $ input
  when
    (len < 2)
    (error $ "convolveGaussian2D: rank = " L.++ show len L.++ " is less than 2.")
  inputF <-
    fmap (fromUnboxed (extent input) . VG.convert) .
    dftExecute
      plan
      (DFTPlanID
         DFT1DG
         (L.reverse . listOfShape . extent $ input)
         [len - 2, len - 1]) .
    VG.convert . toUnboxed . computeS . delay $
    input
  fmap (fromUnboxed (extent input) . VG.convert) .
    dftExecute
      plan
      (DFTPlanID
         IDFT1DG
         (L.reverse . listOfShape . extent $ input)
         [len - 2, len - 1]) .
    VG.convert . toUnboxed . computeS . R.traverse2 inputF filterF const $ \fInput fFilter idx@(sh :. i :. j) ->
    fInput idx * fFilter (Z :. i :. j)
