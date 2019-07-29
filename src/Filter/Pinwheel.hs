{-# LANGUAGE FlexibleContexts #-}
module Filter.Pinwheel
  ( module Filter.Utils
  , PinwheelParams(..)
  , pinwheel
  , pinwheelFilter
  , convolvePinwheel
  , frequencyDomainMultiply
  ) where

import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List           as L (length, reverse)
import           Data.Vector.Generic as VG (convert)
import           DFT.Plan
import           Filter.Utils
import           Utils.Coordinates

data PinwheelParams = PinwheelParams
  { pinwheelParamsRows         :: !Int
  , pinwheelParamsCols         :: !Int
  , pinwheelParamsAlpha        :: !Double
  , pinwheelParamsRadialPeriod :: !Double
  , pinwheelParamsAngularFreqs :: ![Double]
  , pinwheelParamsRadialFreqs  :: ![Double]
  }

{-# INLINE pinwheel #-}
pinwheel :: Double -> Double -> Double -> Double -> Int -> Int -> Complex Double
pinwheel af rf maxR alpha x y
  | r == 0 && af == 0 && rf == 0 = 1
  | r == 0 = 0
  -- | r < 4 = 0
  | otherwise =
    -- (r :+ 0) ** (alpha :+ 0) *
    -- exp (0 :+ rf * ((log r) / (log maxR) * 2 * pi + 0.5 * pi)) *
    -- exp (0 :+ af * theta)
    (r :+ 0) ** (alpha :+ rf * 2 * pi / log maxR) * exp (0 :+ af * theta)
  where
    r = sqrt . fromIntegral $ x ^ (2 :: Int) + y ^ (2 :: Int)
    theta = angleFunctionRad (fromIntegral x) (fromIntegral y)  


{-# INLINE pinwheelPI #-}
pinwheelPI :: Double -> Double -> Double -> Double -> Int -> Int -> Complex Double
pinwheelPI af rf maxR alpha x y
  | r == 0 = 0
  | otherwise =
    (((r) :+ 0) ** (alpha :+ (rf * 2 * pi / (log maxR)))) *
    exp (0 :+ ((af) * theta)) 
  where
    r = (sqrt . fromIntegral $ x ^ (2 :: Int) + y ^ (2 :: Int))
    theta = pi  + angleFunctionRad (fromIntegral x) (fromIntegral y)  

pinwheelFilter ::
     DFTPlan
  -> PinwheelParams
  -> IO (DFTPlan, Array U DIM4 (Complex Double), Array U DIM4 (Complex Double))
pinwheelFilter plan (PinwheelParams rows cols alpha radialPeriod angularFreqs radialFreqs) = do
  let numAngularFreq = L.length angularFreqs
      numRadialFreq = L.length radialFreqs
      filterVec =
        VG.convert . toUnboxed . computeS . makeFilter2D $
        R.traverse2
          (fromListUnboxed (Z :. L.length angularFreqs) angularFreqs)
          (fromListUnboxed (Z :. L.length radialFreqs) radialFreqs)
          (\(Z :. numAngularFreq') (Z :. numRadialFreq') ->
             (Z :. numAngularFreq' :. numRadialFreq' :. cols :. rows))
          (\f1 f2 (Z :. af :. rf :. i :. j) ->
             pinwheel
               (f1 (Z :. af))
               (f2 (Z :. rf))
               radialPeriod
               alpha
               (i - center cols)
               (j - center rows))
      filterPIVec =
        VG.convert . toUnboxed . computeS . makeFilter2D $
        R.traverse2
          (fromListUnboxed (Z :. L.length angularFreqs) angularFreqs)
          (fromListUnboxed (Z :. L.length radialFreqs) radialFreqs)
          (\(Z :. numAngularFreq') (Z :. numRadialFreq') ->
             (Z :. numAngularFreq' :. numRadialFreq' :. cols :. rows))
          (\f1 f2 (Z :. af :. rf :. i :. j) ->
             pinwheelPI
               (f1 (Z :. af))
               (f2 (Z :. rf))
               radialPeriod
               alpha
               (i - div cols 2)
               (j - div rows 2))
      forwardPlanID =
        DFTPlanID DFT1DG [numAngularFreq, numRadialFreq, cols, rows] [2, 3]
      backwardPlanID =
        DFTPlanID IDFT1DG [numAngularFreq, numRadialFreq, cols, rows] [2, 3]
  lock <- getFFTWLock
  (plan1, filterVecF) <-
    dft1dGPlan
      lock
      plan
      [numAngularFreq, numRadialFreq, cols, rows]
      [2, 3]
      filterVec
  (newPlan, _) <-
    idft1dGPlan
      lock
      plan1
      [numAngularFreq, numRadialFreq, cols, rows]
      [2, 3]
      filterVecF
  filterPIVecF <-
    dftExecute
      newPlan
      (DFTPlanID DFT1DG [numAngularFreq, numRadialFreq, cols, rows] [2, 3])
      filterPIVec
  return
    ( newPlan
    , fromUnboxed (Z :. numAngularFreq :. numRadialFreq :. cols :. rows) .
      VG.convert $
      filterVecF
    , fromUnboxed (Z :. numAngularFreq :. numRadialFreq :. cols :. rows) .
      VG.convert $
      filterPIVecF)

{-# INLINE frequencyDomainMultiply #-}
frequencyDomainMultiply ::
     (Source r (Complex Double))
  => Array U DIM4 (Complex Double)
  -> Array r DIM2 (Complex Double)
  -> Array D DIM4 (Complex Double)
frequencyDomainMultiply filterF input =
  R.traverse2
    filterF
    input
    const
    (\f1 f2 idx@(Z :. af :. rf :. i :. j) -> f1 idx * f2 (Z :. i :. j))

{-# INLINE convolvePinwheel #-}
convolvePinwheel ::
     (Source r (Complex Double))
  => DFTPlan
  -> Array U DIM4 (Complex Double)
  -> Array r DIM2 (Complex Double)
  -> IO (Array U DIM4 (Complex Double))
convolvePinwheel plan filterF input =
  fmap (fromUnboxed (extent filterF) . VG.convert) .
  dftExecute
    plan
    (DFTPlanID IDFT1DG (L.reverse . listOfShape . extent $ filterF) [2, 3]) .
  VG.convert . toUnboxed . computeUnboxedS $
  frequencyDomainMultiply filterF input
