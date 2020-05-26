{-# LANGUAGE BangPatterns #-}
module FourierSeries where

import           Control.Monad              as M
import           Data.Array.Repa            as R
import           Data.Complex
import           Data.List                  as L
import           Data.Vector.Storable       as VS
import           Data.Vector.Unboxed        as VU
import           FokkerPlanck.FourierSeries
import           Foreign.CUDA.BLAS          as BLAS
import           Foreign.CUDA.Driver        as CUDA
import           Image.IO
import           Numeric.LinearAlgebra      as NL hiding (All)
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.BLAS
import           Utils.Distribution
import           Utils.Parallel
import           Utils.SimpsonRule
import           Utils.Time

main = do
  args@(gpuIDStr:sizeStr:deltaStr:periodStr:freqStr:angularFreqStr:radialFreqStr:stdStr:_) <-
    getArgs
  printCurrentTime . show $ args
  let gpuID = read gpuIDStr :: Int
      numPoints = read sizeStr :: Int
      delta = read deltaStr :: Double
      period = read periodStr :: Double
      freq = read freqStr :: Int
      numFreq = 2 * freq + 1
      angularFreq = read angularFreqStr :: Int
      radialFreq = read radialFreqStr :: Int
      std = read stdStr :: Double
      folderPath = "output/test/FourierSeries"
      len' = round $ (fromIntegral numPoints / delta :: Double)
      len =
        if odd len'
          then len'
          else len' + 1
      center = div len 2
      !numAngularFreq = 2 * angularFreq + 1
      !numRadialFreq = 2 * radialFreq + 1
      !numPinwheels = numAngularFreq * numRadialFreq
      pinwheel =
        fromFunction (Z :. len :. len :. numAngularFreq :. numRadialFreq) $ \(Z :. i :. j :. af :. rf) ->
          let x = delta * (fromIntegral $ i - center)
              y = delta * (fromIntegral $ j - center)
          in fourierMellin 0.5 (af - angularFreq) (rf - radialFreq) (x, y)
      rectangularHarmonics' =
        fromFunction (Z :. numFreq :. numFreq :. len :. len) $ \(Z :. xFreq :. yFreq :. i :. j) ->
          let !x = fromIntegral $ i - center
              !y = fromIntegral $ j - center
          in cis
               ((-2) * pi * delta *
                ((fromIntegral $ xFreq - freq) * x +
                 (fromIntegral $ yFreq - freq) * y) /
                period)
      simpsonWeights =
        computeUnboxedS $ computeWeightArrFromListOfShape [len, len] :: R.Array U DIM2 (Complex Double)
      inverseR2Harmonics' =
        fromFunction (Z :. numPoints :. numPoints :. numFreq :. numFreq) $ \(Z :. i :. j :. xFreq :. yFreq) ->
          cis
            (2 * pi *
             fromIntegral
               ((xFreq - freq) * (i - (div numPoints 2)) +
                (yFreq - freq) * (j - (div numPoints 2))) /
             period)
      !weight = (delta / 3) ^ 2 :+ 0
  createDirectoryIfMissing True folderPath
  createDirectoryIfMissing True (folderPath </> "Recons")
  createDirectoryIfMissing True (folderPath </> "Coefficients")
  print len
  !weightedPinwheels <-
    computeUnboxedP . R.traverse2 pinwheel simpsonWeights const $ \fP fS idx@(Z :. a :. b :. _ :. _) ->
      fP idx * fS (Z :. a :. b)
  !rectangularHarmonics <- computeUnboxedP rectangularHarmonics'
  !inverseR2Harmonics <- computeUnboxedP inverseR2Harmonics'
  -- printCurrentTime "CPU starts .."
  -- let !gaussianCoefficients =
  --       computeUnboxedS . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
  --         weight *
  --         (gaussian2D
  --            (fromIntegral $ xFreq - freq)
  --            (fromIntegral $ yFreq - freq)
  --            std :+
  --          0)
  -- M.mapM_
  --   (\(af, rf) -> do
  --      let coefficients =
  --            sumS .
  --            sumS .
  --            R.traverse2
  --              rectangularHarmonics
  --              (R.slice weightedPinwheels (Z :. All :. All :. af :. rf))
  --              const $ \fR fP idx@(Z :. _ :. _ :. i :. j) ->
  --              fR idx * fP (Z :. i :. j)
  --          recon =
  --            sumS .
  --            sumS .
  --            R.traverse2
  --              inverseR2Harmonics
  --              (coefficients *^ gaussianCoefficients)
  --              const $ \fR fC idx@(Z :. _ :. _ :. i :. j) ->
  --              fR idx * fC (Z :. i :. j)
  --      plotImageRepaComplex
  --        (folderPath </> "Recon" </> printf "Recon%03d.png" (af * numRadialFreq + rf + 1)) .
  --        ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
  --        recon)
  --   [(af, rf) | af <- [0 .. numAngularFreq - 1], rf <- [0 .. numRadialFreq - 1]]
  -- printCurrentTime "Done"
  -- coefficients <-
  --   (sumP . R.traverse2 rectangularHarmonics weightedPinwheel const $ \fR fP idx@(Z :. _ :. _ :. i :. j) ->
  --      fR idx * fP (Z :. i :. j)) >>=
  --   sumP
  --     gaussianCoefficients =
  --       computeUnboxedS . R.traverse coefficients id $ \f idx@(Z :. xFreq :. yFreq) ->
  --         weight * (f idx) *
  --         (gaussian2D
  --            (fromIntegral $ xFreq - freq)
  --            (fromIntegral $ yFreq - freq)
  --            std :+
  --          0)
  -- recon <-
  --   (sumP . R.traverse2 inverseR2Harmonics gaussianCoefficients const $ \fR fC idx@(Z :. _ :. _ :. i :. j) ->
  --      fR idx * fC (Z :. i :. j)) >>=
  --   sumP
  -- -- Repa
  -- let coefficients =
  --       sumS . sumS . R.traverse2 rectangularHarmonics weightedPinwheel const $ \fR fP idx@(Z :. _ :. _ :. i :. j) ->
  --         fR idx * fP (Z :. i :. j)
  --     gaussianCoefficients =
  --       computeUnboxedS . R.traverse coefficients id $ \f idx@(Z :. xFreq :. yFreq) ->
  --         weight * (f idx) *
  --         (gaussian2D
  --            (fromIntegral $ xFreq - freq)
  --            (fromIntegral $ yFreq - freq)
  --            std :+
  --          0)
  --     recon =
  --       sumS . sumS . R.traverse2 inverseR2Harmonics gaussianCoefficients const $ \fR fC idx@(Z :. _ :. _ :. i :. j) ->
  --         fR idx * fC (Z :. i :. j)
  -- -- HMatrix
  -- let rectangularHarmonicsMat =
  --       ((numFreq ^ 2) >< (len ^ 2)) . R.toList $ rectangularHarmonics
  --     inverseR2HarmonicsMat =
  --       ((numPoints ^ 2) >< (numFreq ^ 2)) . R.toList $ inverseR2Harmonics
  --     weightedPinwheelVec = NL.fromList . R.toList $ weightedPinwheel
  --     gaussianCoefficientsVec =
  --       NL.fromList . R.toList . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
  --         weight *
  --         (gaussian2D
  --            (fromIntegral $ xFreq - freq)
  --            (fromIntegral $ yFreq - freq)
  --            std :+
  --          0)
  --     coefficientsVec = rectangularHarmonicsMat #> weightedPinwheelVec
  --     recon =
  --       fromListUnboxed (Z :. numPoints :. numPoints) . NL.toList $
  --       inverseR2HarmonicsMat #> (coefficientsVec * gaussianCoefficientsVec)
  -- --BLAS mxv
  -- let rectangularHarmonicsMat = VU.convert . toUnboxed $ rectangularHarmonics
  --     inverseR2HarmonicsMat = VU.convert . toUnboxed $ inverseR2Harmonics
  --     weightedPinwheelVec = VU.convert . toUnboxed $ weightedPinwheel
  --     gaussianCoefficientsVec =
  --       VU.convert .
  --       toUnboxed . computeS . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
  --         weight *
  --         (gaussian2D
  --            (fromIntegral $ xFreq - freq)
  --            (fromIntegral $ yFreq - freq)
  --            std :+
  --          0)
  -- coefficientsVec <-
  --   mxv (numFreq ^ 2) (len ^ 2) rectangularHarmonicsMat weightedPinwheelVec
  -- reconVec <-
  --   mxv
  --     (numPoints ^ 2)
  --     (numFreq ^ 2)
  --     inverseR2HarmonicsMat
  --     (VS.zipWith (*) coefficientsVec gaussianCoefficientsVec)
  -- let recon = fromUnboxed (Z :. numPoints :. numPoints) . VS.convert $ reconVec
  --BLAS mxm
  -- let !rectangularHarmonicsMat = VU.convert . toUnboxed $ rectangularHarmonics
  --     !inverseR2HarmonicsMat = VU.convert . toUnboxed $ inverseR2Harmonics
  --     !weightedPinwheelVec = VU.convert . toUnboxed $ weightedPinwheels
  --     !gaussianCoefficients =
  --       computeUnboxedS . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
  --         weight *
  --         (gaussian2D
  --            (fromIntegral $ xFreq - freq)
  --            (fromIntegral $ yFreq - freq)
  --            std :+
  --          0)
  -- printCurrentTime "Start.."
  -- coefficientsVec <-
  --   gemmBLAS
  --     (numFreq ^ 2)
  --     numPinwheels
  --     (len ^ 2)
  --     rectangularHarmonicsMat
  --     weightedPinwheelVec
  -- printCurrentTime "coefficients done"
  -- let coefficientsArr =
  --       fromUnboxed (Z :. numFreq :. numFreq :. numAngularFreq :. numRadialFreq) .
  --       VS.convert $
  --       coefficientsVec
  --     weightedCoefficientsArr =
  --       computeS .
  --       R.traverse2
  --         coefficientsArr
  --         gaussianCoefficients
  --         const $ \f1 f2 idx@(Z :. i :. j :. _ :. _) ->
  --         f1 idx * f2 (Z :. i :. j)
  --     weightedCoefficients = VS.convert . toUnboxed $ weightedCoefficientsArr
  -- M.mapM_
  --   (\(af, rf) ->
  --      plotImageRepaComplex
  --        (folderPath </> "Coefficients" </>
  --         printf
  --           "Coef(%d,%d).png"
  --           (af - div numAngularFreq 2)
  --           (rf - div numRadialFreq 2)) .
  --      ImageRepa 8 .
  --      computeS .
  --      extend (Z :. (1 :: Int) :. All :. All) . R.slice coefficientsArr $
  --      (Z :. All :. All :. af :. rf))
  --   [(af, rf) | af <- [0 .. numAngularFreq - 1], rf <- [0 .. numRadialFreq - 1]]
  -- printCurrentTime "Done"
  -- reconVec <-
  --   gemmBLAS
  --     (numPoints ^ 2)
  --     numPinwheels
  --     (numFreq ^ 2)
  --     inverseR2HarmonicsMat
  --     weightedCoefficients
  -- let !recon =
  --       fromUnboxed
  --         (Z :. numPoints :. numPoints :. numAngularFreq :. numRadialFreq) .
  --       VS.convert $
  --       reconVec
  -- M.mapM_
  --   (\(af, rf) ->
  --      plotImageRepaComplex
  --        (folderPath </> "Recons" </>
  --         printf
  --           "Recon(%d,%d).png"
  --           (af - div numAngularFreq 2)
  --           (rf - div numRadialFreq 2)) .
  --      ImageRepa 8 .
  --      computeS . extend (Z :. (1 :: Int) :. All :. All) . R.slice recon $
  --      (Z :. All :. All :. af :. rf))
  --   [(af, rf) | af <- [0 .. numAngularFreq - 1], rf <- [0 .. numRadialFreq - 1]]
  -- printCurrentTime "Done"
  --BLAS mxmGPU
  let !rectangularHarmonicsMat = VU.convert . toUnboxed $ rectangularHarmonics
      !inverseR2HarmonicsMat = VU.convert . toUnboxed $ inverseR2Harmonics
      !weightedPinwheelVec = VU.convert . toUnboxed $ weightedPinwheels
      !gaussianCoefficients =
        computeUnboxedS . R.fromFunction (Z :. numFreq :. numFreq) $ \(Z :. xFreq :. yFreq) ->
          weight *
          (gaussian2D
             (fromIntegral $ xFreq - freq)
             (fromIntegral $ yFreq - freq)
             std :+
           0)
  initialise []
  dev <- device gpuID
  ctx <- CUDA.create dev []
  handle <- BLAS.create
  printCurrentTime "Start.."
  coefficientsVec <-
    gemmCuBLAS11
      handle
      (numFreq ^ 2)
      numPinwheels
      (len ^ 2)
      rectangularHarmonicsMat
      weightedPinwheelVec
  let coefficientsArr =
        fromUnboxed (Z :. numFreq :. numFreq :. numAngularFreq :. numRadialFreq) .
        VS.convert $
        coefficientsVec
      weightedCoefficients =
        VU.convert .
        toUnboxed .
        computeS . R.traverse2 coefficientsArr gaussianCoefficients const $ \f1 f2 idx@(Z :. i :. j :. _ :. _) ->
          f1 idx * f2 (Z :. i :. j)
  printCurrentTime "coefficients done"
  reconVec <-
    gemmCuBLAS11
      handle
      (numPoints ^ 2)
      numPinwheels
      (numFreq ^ 2)
      inverseR2HarmonicsMat
      weightedCoefficients
  let !recon =
        fromUnboxed (Z :. numPoints :. numPoints :. numPinwheels) . VS.convert $
        reconVec
  M.mapM_
    (\i ->
       plotImageRepaComplex
         (folderPath </> "Recons" </> printf "Recon%03d.png" (i + 1)) .
       ImageRepa 8 .
       computeS . extend (Z :. (1 :: Int) :. All :. All) . R.slice recon $
       (Z :. All :. All :. i))
    [0 .. numPinwheels - 1]
  printCurrentTime "Done"
  M.mapM_
    (\(af, rf) ->
       plotImageRepaComplex
         (folderPath </> "Coefficients" </>
          printf
            "Coef(%d,%d).png"
            (af - div numAngularFreq 2)
            (rf - div numRadialFreq 2)) .
       ImageRepa 8 .
       computeS .
       extend (Z :. (1 :: Int) :. All :. All) . R.slice coefficientsArr $
       (Z :. All :. All :. af :. rf))
    [(af, rf) | af <- [0 .. numAngularFreq - 1], rf <- [0 .. numRadialFreq - 1]]
  -- plotImageRepaComplex (folderPath </> "Grating.png") .
  --   ImageRepa 8 .
  --   computeS .
  --   extend (Z :. (1 :: Int) :. All :. All) . R.slice rectangularHarmonics $
  --   (Z :. (0 :: Int) :. freq :. All :. All)
  -- plotImageRepaComplex (folderPath </> "Pinwheel.png") .
  --   ImageRepa 8 .
  --   computeS . fromFunction (Z :. (1 :: Int) :. numPoints :. numPoints) $ \(Z :. _ :. i :. j) ->
  --   fourierMellin
  --     angularFreq
  --     radialFreq
  --     ( (fromIntegral $ i - (div numPoints 2))
  --     , (fromIntegral $ j - (div numPoints 2)))
