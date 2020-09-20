{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Strict #-}
module FokkerPlanck.GreensFunction where

import           Control.DeepSeq
import           Control.Monad                  as M
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.IArray              as IA
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           FokkerPlanck.Analytic
import           FokkerPlanck.FourierSeriesGPU
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.FourierSeries2D
import           Image.IO
import           Pinwheel.FourierSeries2D
import           Pinwheel.Transform
import           Sparse.Vector
import           STC.Utils
import           System.FilePath
import           Utils.Array
import           Utils.List
import           Utils.Parallel
import           Utils.SimpsonRule
import           Utils.Time
import Text.Printf
import           Utils.Distribution

sampleCartesian ::
     FilePath
  -> FilePath
  -> [PTX]
  -> Int
  -> Double
  -> Double
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> Int
  -> Int
  -> Int -> Double
  -> IO (Histogram (Complex Double))
sampleCartesian filePath folderPath ptxs numPoints period delta oris gamma thetaSigma tau threshold r2Sigma maxPhiFreqs maxRhoFreqs maxThetaFreqs maxRFreqs stdR2 = do
  let center = div numPoints 2
      deltaTheta = 2 * pi / fromIntegral oris
      origin = R2S1RP 0 0 0 gamma
      idx =
        L.filter
          (\(R2S1RP c r _ g) -> (sqrt (c ^ 2 + r ^ 2) > 0) && g > 0)
          [ R2S1RP
            ((fromIntegral $ c - center) * delta)
            ((fromIntegral $ r - center) * delta)
            0
            gamma
          | c <- [0 .. numPoints - 1]
          , r <- [0 .. numPoints - 1]
          ]
      logGamma = 0 -- log gamma
      -- simpsonWeights = weightsSimpsonRule oris
      xs =
        L.concat -- .
        -- L.zipWith
        --   (\w -> L.map (\(a, b, c, d, v) -> (a, b, c, d, v * w)))
        --   simpsonWeights
         $
        parMap
          rdeepseq
          (\o ->
             let ori = fromIntegral o * deltaTheta
              in L.filter (\(_, _, _, _, v) -> v > threshold) .
                 L.map
                   (\point@(R2S1RP x y _ g) ->
                      let v =
                            computePji thetaSigma tau origin (R2S1RP x y ori g)
                          phi = atan2 y x
                          rho = sqrt $ (x ^ 2 + y ^ 2)
                       in ( phi
                          , log rho
                          , ori
                          , logGamma
                          , v * (1 - gaussian2DPolar rho stdR2)
                          )) $
                 idx)
          [0 .. oris - 1]
  -- arr <-
  --   computeUnboxedP . fromFunction (Z :. oris :. numPoints :. numPoints) $ \(Z :. o :. i :. j) ->
  --     let x = (fromIntegral $ i - center) * delta
  --         y = (fromIntegral $ j - center) * delta
  --         phi = atan2 y x
  --         rho = sqrt $ (x ^ 2 + y ^ 2)
  --         theta = fromIntegral o * deltaTheta
  --         v = computePji thetaSigma tau origin (R2S1RP x y theta gamma)
  --     in if rho <= 1.5
  --          then (1, 1, 1, 1, 0)
  --          else (phi, log rho, theta, logGamma, v)
  -- plotImageRepa (folderPath </> "GreensFunction.png") .
  --   ImageRepa 8 .
  --   computeS .
  --   -- reduceContrast 20 .
  --   extend (Z :. (1 :: Int) :. All :. All) .
  --   sumS . rotate3D . R.map (\(_, _, _, _, v) -> v) $
  --   arr
      -- simpsonWeights =
      --   computeWeightArrFromListOfShape [numPoints, numPoints, oris]
      -- xs =
      --   L.filter (\(_, _, _, _, v) -> v > threshold) .
      --   R.toList -- .
      --   -- R.zipWith (\w (a, b, c, d, v) -> (a, b, c, d, w * v)) simpsonWeights
      --   $
      --   arr
  let
  printCurrentTime $
    printf
      "Sparsity: %f%%\n"
      (100 * (fromIntegral . L.length $ xs) /
       (fromIntegral $ oris * numPoints ^ 2) :: Double)
  let phiFreqs = L.map fromIntegral [-maxPhiFreqs .. maxPhiFreqs]
      rhoFreqs = L.map fromIntegral [-maxRhoFreqs .. maxRhoFreqs]
      thetaFreqs = L.map fromIntegral [-maxThetaFreqs .. maxThetaFreqs]
      rFreqs = L.map fromIntegral [-maxRFreqs .. maxRFreqs]
      hist =
        L.foldl1' (addHistogram) .
        parZipWith
          rdeepseq
          (\ptx ys ->
             computeFourierCoefficientsGPU'
               r2Sigma
               period
               phiFreqs
               rhoFreqs
               thetaFreqs
               rFreqs
               ptx
               ys)
          ptxs .
        divideListN (L.length ptxs) $
        xs
  encodeFile filePath hist
  return hist


-- sampleLogpolar ::
--      FilePath
--   -> [PTX]
--   -> Int
--   -> Int
--   -> Int
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Double
--   -> Int
--   -> Int
--   -> Int
--   -> Int
--   -> IO (Histogram (Complex Double))
-- sampleLogpolar filePath ptxs phis rhos oris maxScale gamma sigma tau threshold r2Sigma maxPhiFreqs maxRhoFreqs maxThetaFreqs maxRFreqs = do
--   let deltaPhi = 2 * pi / fromIntegral phis
--       deltaTheta = 2 * pi / fromIntegral oris
--       logMaxScale = log maxScale
--       deltaLogRho = logMaxScale / fromIntegral rhos
--       logGamma = 0 -- log gamma
--       origin = R2S1RP 0 0 0 gamma
--       -- xs =
--       --   parMap
--       --     rdeepseq
--       --     (\(phi, logRho, theta) ->
--       --        let x = (exp logRho) * cos phi
--       --            y = (exp logRho) * sin phi
--       --            v = computePji sigma tau origin (R2S1RP x y theta gamma)
--       --        in (phi, logRho, theta, 0, v))
--       --     [ ( deltaTheta * fromIntegral iPhi
--       --       , deltaLogRho * fromIntegral iLogRho - logMaxScale
--       --       , deltaTheta * fromIntegral iTheta)
--       --     | iPhi <- [0 .. oris - 1]
--       --     , iLogRho <- [0 .. rhos - 1]
--       --     , iTheta <- [0 .. oris - 1]
--       --     ]
--   arr <-
--     computeUnboxedP . R.fromFunction (Z :. oris :. rhos :. phis) $ \(Z :. iTheta :. iLogRho :. iPhi) ->
--       let theta = deltaTheta * fromIntegral iTheta
--           phi = deltaPhi * fromIntegral iPhi
--           logRho = deltaLogRho * fromIntegral iLogRho -- logMaxScale
--           x = (exp logRho) * cos phi
--           y = (exp logRho) * sin phi
--           v = computePji sigma tau origin (R2S1RP x y theta gamma)
--       in (phi, logRho, theta, logGamma, v)
--   let simpsonWeights = computeWeightArrFromListOfShape [phis, rhos, oris]
--       xs =
--         L.filter (\(_, _, _, _, v) -> v > threshold) .
--         R.toList .
--         R.zipWith (\w (a, b, c, d, v) -> (a, b, c, d, w * v)) simpsonWeights $
--         arr
--   printCurrentTime $
--     printf
--       "Sparsity: %f%%\n"
--       (100 * (fromIntegral . L.length $ xs) / (fromIntegral $ oris * phis * rhos) :: Double)
--   let phiFreqs = L.map fromIntegral [-maxPhiFreqs .. maxPhiFreqs]
--       rhoFreqs = L.map fromIntegral [-maxRhoFreqs .. maxRhoFreqs]
--       thetaFreqs = L.map fromIntegral [-maxThetaFreqs .. maxThetaFreqs]
--       rFreqs = L.map fromIntegral [-maxRFreqs .. maxRFreqs]
--       hist =
--         L.foldl1' (addHistogram) .
--         parZipWith
--           rdeepseq
--           (\ptx ys ->
--              computeFourierCoefficientsGPU'
--                r2Sigma
--                phiFreqs
--                rhoFreqs
--                thetaFreqs
--                rFreqs
--                ptx
--                ys)
--           ptxs .
--         divideListN (L.length ptxs) $
--         xs
--   encodeFile filePath hist
--   return hist

{-# INLINE sampleR2S1 #-}
sampleR2S1 ::
     Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Int]
  -> IO (R.Array U DIM3 Double)
sampleR2S1 rows cols delta deltaTheta gamma sigma tau oris = do
  let centerRow = div rows 2
      centerCol = div cols 2
      origin = R2S1RP 0 0 0 gamma
  computeUnboxedP .
    R.traverse
      (fromListUnboxed (Z :. (L.length oris)) oris)
      (\(Z :. o) -> (Z :. cols :. rows :. o)) $ \f (Z :. c :. r :. o) ->
    if (c == centerCol) && (r == centerRow)
      then 0
      else computePji
             sigma
             tau
             origin
             (R2S1RP
                (delta * fromIntegral (c - centerCol))
                (delta * fromIntegral (r - centerRow))
                (fromIntegral (f (Z :. o)) * deltaTheta)
                gamma)

{-# INLINE sampleR2S1' #-}
sampleR2S1' ::
     Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> R.Array U DIM2 (Int, Double)
sampleR2S1' rows cols delta gamma sigma tau ori =
  let centerRow = div rows 2
      centerCol = div cols 2
      origin = R2S1RP 0 0 0 gamma
  in computeUnboxedS . R.fromFunction (Z :. cols :. rows) $ \(Z :. c :. r) ->
       let x = c - centerCol
           y = r - centerRow
       in if x == 0 && y == 0
            then (0, 0)
            else ( getVectorIndex2D rows c r
                 , computePji
                     sigma
                     tau
                     origin
                     (R2S1RP
                        (delta * fromIntegral (c - centerCol))
                        (delta * fromIntegral (r - centerRow))
                        ori
                        gamma))

computeFourierCoefficients ::
     FilePath
  -> [Int]
  -> [PTX]
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> Double
  -> Double
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> IO (Histogram (Complex Double))
computeFourierCoefficients histFilePath deviceIDs ptxs numPoints oris delta deltaFreq gamma sigma tau numR2Freqs period std numBatchR2Freqs numBatchOri maxPhiFreq maxRhoFreq maxThetaFreq maxRFreq = do
  when
    (even oris)
    (error "computeFourierCoefficients: the number of orientations is not odd.")
  let deltaTheta = 2 * pi / fromIntegral oris
      simpsonWeights =
        computeUnboxedS $
        (computeWeightArrFromListOfShape [numPoints, numPoints] :: R.Array D DIM2 (Complex Double))
      oriIdxs = divideListN numBatchOri [0 .. oris - 1]
  xs <-
    M.mapM
      (sampleR2S1 numPoints numPoints delta deltaTheta gamma sigma tau)
      oriIdxs
  let weightedVecs =
        parMap
          rdeepseq
          (\arr ->
             let (Z :. n :. m :. c) = extent arr
             in CuMat (n * m) c .
                CuVecHost .
                VU.convert .
                toUnboxed . computeS . R.traverse2 arr simpsonWeights const $ \f fW idx@(Z :. i :. j :. _) ->
                  (f idx :+ 0) * fW (Z :. i :. j))
          xs
  coefficients <-
    computeFourierCoefficientsR2Stream
      deviceIDs
      ptxs
      numR2Freqs
      numPoints
      period
      delta
      deltaFreq
      numBatchR2Freqs
      weightedVecs
  print "coefficients done"
  let centerR2Freq = div numR2Freqs 2
      centerPoints = div numPoints 2
      logGamma = log gamma
      gaussianWeightedCoef =
        computeUnboxedS $ applyGaussian deltaFreq std coefficients
      -- gaussianWeightedCoef = coefficients
      simpsonWeightsOri =
        computeUnboxedS $
        (computeWeightArrFromListOfShape [oris] :: R.Array D DIM1 (Complex Double))
      weightedCoef =
        R.traverse2 gaussianWeightedCoef simpsonWeightsOri const $ \fG fS idx@(Z :. o :. _ :. _) ->
          fG idx * fS (Z :. o)
      folderPath = "output/test/STCPinwheel"
  -- print . extent $ gaussianWeightedCoef
  -- inverseR2Harmonics <-
  --   createInverseHarmonicMatriesGPU ptxs 1 numPoints numR2Freqs period delta
  -- sourceR2 <-
  --   computeFourierSeriesR2
  --     deviceIDs
  --     numR2Freqs
  --     numPoints
  --     period
  --     inverseR2Harmonics
  --     [ CuMat oris (numR2Freqs ^ 2) . CuVecHost . VU.convert . toUnboxed $
  --       gaussianWeightedCoef
  --     ]
  plotImageRepa (folderPath </> "Green.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (x) ** 2) . L.head $
    xs
  plotImageRepa (folderPath </> "Coefficients.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freqs :. numR2Freqs) .
    VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
    gaussianWeightedCoef
  -- plotImageRepa (folderPath </> "GreenRecon.png") .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
  --   VU.map sqrt . toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
  --   sourceR2
  ys <-
    (computeUnboxedP . R.traverse weightedCoef id $ \f idx@(Z :. o :. i :. j) ->
       if i == centerR2Freq && j == centerR2Freq
         then (0, 0, 0, 0, 0)
         else let x = deltaFreq * (fromIntegral $ i - centerR2Freq)
                  y = deltaFreq * (fromIntegral $ j - centerR2Freq)
                  phi = atan2 y x
                  logRho = log . sqrt $ (x ^ 2 + y ^ 2)
              in (phi, logRho, fromIntegral o * deltaTheta, logGamma, f idx))
  -- ys <-
  --   (computeUnboxedP .
  --    R.traverse
  --      (R.traverse2 (L.head xs) simpsonWeightsOri const $ \fG fS idx@(Z :. _ :. _ :. o) ->
  --         (fG idx :+ 0) * fS (Z :. o))
  --      id $ \f idx@(Z :. i :. j :. o) ->
  --      if i == centerPoints && j == centerPoints
  --        then (0, 0, 0, 0, 0)
  --        else let x = delta * (fromIntegral $ i - centerPoints)
  --                 y = delta * (fromIntegral $ j - centerPoints)
  --                 phi = atan2 x y
  --                 logRho = log . sqrt $ (x ^ 2 + y ^ 2)
  --             in ( phi
  --                , logRho
  --                , fromIntegral o * deltaTheta
  --                , logGamma
  --                , f idx))
  let hist =
        L.foldl1' (addHistogram) .
        parZipWith
          rdeepseq
          (\ptx ->
             computePinwheelCoefficients
               maxPhiFreq
               maxRhoFreq
               maxThetaFreq
               maxRFreq
               ptx)
          ptxs .
        divideListN (L.length ptxs) . R.toList $
        ys
  encodeFile histFilePath hist
  return hist



computePinwheelTransformCoefficients ::
     Int
  -> FilePath
  -> [Int]
  -> [PTX]
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> Double
  -> Double
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> IO (Histogram (Complex Double))
computePinwheelTransformCoefficients numThread histFilePath deviceIDs ptxs numPoints oris delta threshold gamma sigma tau numR2Freqs periodR2 std numBatchR2 numBatchPinwheelFreqs maxPhiFreq maxRhoFreq maxThetaFreq maxRFreq batchSize s = do
  printCurrentTime "compute pinwheelArr"
  pinwheelArr <-
    pinwheelFourierSeries
      deviceIDs
      ptxs
      numPoints -- numR2Freqs
      numPoints
      delta
      periodR2
      maxPhiFreq
      maxRhoFreq
      maxThetaFreq
      maxRFreq
      s
      numBatchR2
      numBatchPinwheelFreqs :: IO (IA.Array (Int, Int) (VU.Vector (Complex Double)))
  printCurrentTime "pinwheelArr done"
  let folderPath = "output/test/STCPinwheel"
  let deltaTheta = 2 * pi / (fromIntegral oris)
      greensFunc =
        parMap
          rdeepseq
          (\o ->
             let theta = fromIntegral o * deltaTheta
                 vec =
                   toUnboxed $
                   sampleR2S1' numPoints numPoints delta gamma sigma tau theta
             in PinwheelTransformData (log gamma) theta .
                uncurry SparseVector .
                VU.unzip .
                VU.map (\(i, v) -> (i, v :+ 0)) .
                VU.filter (\(_, v) -> v > threshold) $
                vec)
          [0 .. oris - 1]
      len =
        L.foldl'
          (\s (PinwheelTransformData _ _ (SparseVector vec _)) ->
             s + VU.length vec)
          0
          greensFunc
  printCurrentTime $
    printf
      "Sparsity: %f%%\n"
      (100 * (fromIntegral len) / (fromIntegral $ oris * numPoints ^ 2) :: Double)
  let hist =
        pinwheelTransform
          numThread
          maxPhiFreq
          maxRhoFreq
          maxThetaFreq
          maxRFreq
          s
          pinwheelArr
          greensFunc
  encodeFile histFilePath hist
  return hist
