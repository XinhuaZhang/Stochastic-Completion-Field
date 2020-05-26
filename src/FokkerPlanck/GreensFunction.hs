{-# LANGUAGE BangPatterns #-}
module FokkerPlanck.GreensFunction where

import           Control.Monad                  as M
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           FokkerPlanck.Analytic
import           FokkerPlanck.FourierSeriesGPU
import           FokkerPlanck.Histogram
import           Foreign.CUDA.Driver            as CUDA
import           Image.IO
import           STC.Utils
import           System.FilePath
import           Utils.Array
import           Utils.List
import           Utils.Parallel

sampleCartesian ::
     FilePath
  -> FilePath
  -> [Int]
  -> Double
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> IO (Histogram (Complex Double))
sampleCartesian !folderPath !filePath !deviceID !radius !oris !delta !gamma !sigma !tau !phiFreqs !rhoFreqs !thetaFreqs !rFreqs = do
  let !rows = round $ (2 * radius + 1) / delta
      !cols = rows
      !centerRow = fromIntegral $ div rows 2
      !centerCol = fromIntegral $ div cols 2
      !deltaTheta = 2 * pi / fromIntegral oris
      !origin = R2S1RP 0 0 0 gamma
      !idx =
        L.filter
          (\(R2S1RP c r _ g) -> ((c /= 0) || (r /= 0)) && g > 0)
          [ R2S1RP
            ((fromIntegral c - centerCol) * delta)
            ((fromIntegral r - centerRow) * delta)
            0
            gamma
          | c <- [0 .. cols - 1]
          , r <- [0 .. rows - 1]
          ]
      !logGamma = log gamma
      !xs =
        L.concat $
        parMap
          rdeepseq
          (\o ->
             let !ori = fromIntegral o * deltaTheta
             in L.filter (\(_, _, _, _, v) -> v > 1e-8) .
                L.map
                  (\point@(R2S1RP x y _ g) ->
                     let !v = computePji sigma tau origin (R2S1RP x y ori g)
                         !phi = atan2 y x
                         !r = log . sqrt $ (x ^ 2 + y ^ 2)
                     in  (phi, r, ori, logGamma, v)) $
                idx)
          [0 .. oris - 1]
      len = L.length xs
      s = L.sum . L.map (\(_, _, _, _, v) -> v) $ xs
      -- arr = fromListUnboxed (Z :. oris :. cols :. rows) xs
  -- plotImageRepa
  --   (folderPath </> "SourceR2S1.png")
  --   (ImageRepa 8 .
  --    computeS .
  --    R.extend (Z :. (1 :: Int) :. All :. All) .
  --    -- reduceContrast 100 .
  --    R.sumS . rotate3D . R.map (\(_, _, _, _, v) -> v) $
  --    arr)
  print (s / fromIntegral len)
  initialise []
  devs <- M.mapM device deviceID
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let !hist =
        L.foldl1' (addHistogram) .
        parZipWith
          rdeepseq
          (\ptx ys ->
             computeFourierCoefficientsGPU'
               phiFreqs
               rhoFreqs
               thetaFreqs
               rFreqs
               ptx
               ys)
          ptxs .
        divideList (div (L.length xs) (L.length deviceID)) $
        xs
  encodeFile filePath hist
  return hist


sampleLogpolar ::
     FilePath
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> IO (Histogram (Complex Double))
sampleLogpolar !filePath !deviceID !rhos !oris !maxScale !gamma !sigma !tau !phiFreqs !rhoFreqs !thetaFreqs !rFreqs = do
  let !deltaTheta = 2 * pi / fromIntegral oris
      !logMaxScale = log maxScale
      !deltaLogRho = 2 * logMaxScale / fromIntegral rhos
      !origin = R2S1RP 0 0 0 gamma
      !xs =
        parMap
          rdeepseq
          (\(phi, logRho, theta) ->
             let !x = (exp logRho) * cos phi
                 !y = (exp logRho) * sin phi
                 !v = computePji sigma tau origin (R2S1RP x y theta gamma)
             in (phi, logRho, theta, 0, v))
          [ ( deltaTheta * fromIntegral iPhi
            , deltaLogRho * fromIntegral iLogRho - logMaxScale
            , deltaTheta * fromIntegral iTheta)
          | iPhi <- [0 .. oris - 1]
          , iLogRho <- [0 .. rhos - 1]
          , iTheta <- [0 .. oris - 1]
          ]
  initialise []
  dev <- device deviceID
  ctx <- CUDA.create dev []
  ptx <- createTargetFromContext ctx
  let !hist =
        computeFourierCoefficientsGPU'
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          ptx
          xs
  encodeFile filePath hist
  return hist


sampleR2S1 ::
     Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Double
  -> Double
  -> R.Array U DIM3 Double
sampleR2S1 !rows !cols !oris !delta !gamma !sigma !tau =
  let !centerRow = fromIntegral $ div rows 2
      !centerCol = fromIntegral $ div cols 2
      !deltaTheta = 2 * pi / fromIntegral oris
      !origin = R2S1RP 0 0 0 gamma
      !weights =
        [1 / 16, 1 / 8, 1 / 16, 1 / 8, 1 / 4, 1 / 8, 1 / 16, 1 / 8, 1 / 16]
      -- !origins =
      --   L.zipWith
      --     (\w (i, j) -> (R2S1RP (i * delta) (j * delta) 0 gamma, w))
      --     weights
      --     [(i, j) | i <- [-1, 0, 1], j <- [-1, 0, 1]]
      -- !ys =
      --   L.zipWith
      --     (\w (i, j) -> (i * delta, j * delta, w))
      --     weights
      --     [(i, j) | i <- [-1, 0, 1], j <- [-1, 0, 1]]
      idx =
        [ R2S1RP
          ((fromIntegral c - centerCol) * delta)
          ((fromIntegral r - centerRow) * delta)
          0
          gamma
        | c <- [0 .. cols - 1]
        , r <- [0 .. rows - 1]
        ]
      !xs =
        L.concat $
        parMap
          rdeepseq
          (\o ->
             L.map
               (\point@(R2S1RP x y _ g) ->
                  if (x /= 0) || (y /= 0)
                    -- then L.foldl'
                    --        (\s (x', y', w) ->
                    --           s +
                    --           if x + x' == 0 && y + y' == 0
                    --             then 0
                    --             else w *
                    --                  computePji
                    --                    sigma
                    --                    tau
                    --                    origin
                    --                    (R2S1RP
                    --                       (x + x')
                    --                       (y + y')
                    --                       (fromIntegral o * deltaTheta)
                    --                       g))
                    --        0
                    --        ys
                    then computePji
                           sigma
                           tau
                           origin
                           (R2S1RP x y (fromIntegral o * deltaTheta) g)
                    else 0)
               idx)
          [0 .. oris - 1]
  in fromListUnboxed (Z :. oris :. cols :. rows) xs
