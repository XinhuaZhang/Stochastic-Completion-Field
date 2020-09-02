{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
module FokkerPlanck.FourierSeriesGPU where

import           Control.Monad                      as M
import           Control.Monad.Parallel             as MP
import qualified Data.Array.Accelerate              as A
import qualified Data.Array.Accelerate.Data.Complex as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Complex                       as C
import           Data.DList                         as DL
import           Data.List                          as L
import           Data.Vector.Unboxed                as VU
import           FokkerPlanck.Analytic
import           FokkerPlanck.BrownianMotion        (Particle (..),
                                                     generateRandomNumber,
                                                     moveParticle, scalePlus,
                                                     thetaPlus)
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.GPUKernel
import           FokkerPlanck.Histogram
import           GHC.Float
import           Statistics.Distribution
import           Statistics.Distribution.Normal
import           System.Random.MWC
import           Utils.Parallel
import Debug.Trace
import Utils.List

{-# INLINE computeFourierCoefficientsGPU #-}
computeFourierCoefficientsGPU ::
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Acc (A.Vector (Float, Float, Float, Float))
  -> A.Exp Float
  -> A.Exp Float
  -> PTX
  -> [DList Particle]
  -> Histogram (C.Complex Double)
computeFourierCoefficientsGPU !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !freqArr !sigmaExp !logPeriodExp !ptx !xs =
  let particles =
        L.map
          (\(Particle phi rho theta r w) ->
             ( double2Float phi
             , double2Float (log rho)
             , double2Float theta
             , double2Float (log r)
             , double2Float w)) .
        DL.toList . DL.concat $
        xs
      !len = L.length particles
   in Histogram
        [ L.length phiFreqs
        , L.length rhoFreqs
        , L.length thetaFreqs
        , L.length rFreqs
        ]
        len .
      VU.map (\(a :+ b) -> float2Double a :+ float2Double b) .
      VU.fromList .
      A.toList .
      runNWith ptx (gpuKernel'' sigmaExp logPeriodExp freqArr) .
      A.fromList (A.Z A.:. len) $
      particles


{-# INLINE computeFrequencyArray #-}
computeFrequencyArray ::
     (A.Elt a) => [a] -> [a] -> [a] -> [a] -> A.Vector (a, a, a, a)
computeFrequencyArray !phiFreqs !rhoFreqs !thetaFreqs !rFreqs =
  A.fromList
    (A.Z A.:.
     ((L.length rFreqs) * (L.length thetaFreqs) * (L.length rhoFreqs) *
      L.length phiFreqs))
    [ (rFreq', thetaFreq', rhoFreq', phiFreq')
    | rFreq' <- rFreqs
    , thetaFreq' <- thetaFreqs
    , rhoFreq' <- rhoFreqs
    , phiFreq' <- phiFreqs
    ]

-- {-# INLINE computeFourierCoefficientsGPU #-}
-- computeFourierCoefficientsGPU ::
--      Acc (A.Vector (Double, Double, Double, Double))
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> Double
--   -> Exp Double
--   -> Exp Double
--   -> Exp (AC.Complex Double)
--   -> PTX
--   -> [DList Particle]
--   -> Histogram (C.Complex Double)
-- computeFourierCoefficientsGPU !freqArr !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !deltaLogRho !maxScaleExp !halfLogPeriodExp !deltaLogRhoComplexExp !ptx !xs =
--   let !particles = DL.toList . DL.concat $ xs
--       !len = L.length particles
--       !scaleSampledParticles =
--         L.concat .
--         parMap
--           rdeepseq
--           (\particle@(Particle phi rho theta r) ->
--              let !n = Prelude.floor $ (log r + halfLogPeriod) / deltaLogRho
--              in (phi, rho, theta, r) :
--                 [ ( phi
--                   , rho
--                   , theta
--                   , (exp $
--                      (Prelude.fromIntegral i * deltaLogRho - halfLogPeriod)))
--                 | i <- [0 .. n - 1]
--                 ]) $
--         particles
--       !movedParticleArr =
--         compute .
--         afst .
--         A.filter
--           (\particle ->
--              let (_, rho, _, _) =
--                    unlift particle :: ( Exp Double
--                                       , Exp Double
--                                       , Exp Double
--                                       , Exp Double)
--              in rho A./= 0) .
--         A.map (moveParticle maxScaleExp) .
--         use . A.fromList (Z :. (L.length scaleSampledParticles)) $
--         scaleSampledParticles
--   in Histogram
--        [ L.length phiFreqs
--        , L.length rhoFreqs
--        , L.length thetaFreqs
--        , L.length rFreqs
--        ]
--        len .
--      VU.fromList . A.toList $
--      runWith ptx .
--      A.map
--        (\(unlift -> (rFreq, thetaFreq, rhoFreq, phiFreq)) ->
--           sfoldl
--             (\s particle ->
--                s +
--                deltaLogRhoComplexExp *
--                (coefficient
--                   halfLogPeriodExp
--                   rFreq
--                   thetaFreq
--                   rhoFreq
--                   phiFreq
--                   particle))
--             0
--             (constant Z)
--             movedParticleArr) $
--      freqArr


computeFourierCoefficientsGPU' ::
     Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> PTX
  -> [(Double, Double, Double, Double, Double)]
  -> Histogram (Complex Double)
computeFourierCoefficientsGPU' !sigma !period !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !ptx !xs =
  let -- !freqArr =
      --   A.use $
      --   computeFrequencyArray
      --     (L.map double2Float phiFreqs)
      --     (L.map double2Float rhoFreqs)
      --     (L.map double2Float thetaFreqs)
      --     (L.map double2Float rFreqs)
      !freqArr =
        A.use $
        computeFrequencyArray
          ( phiFreqs)
          ( rhoFreqs)
          ( thetaFreqs)
          ( rFreqs)
   in Histogram
        [ L.length phiFreqs
        , L.length rhoFreqs
        , L.length thetaFreqs
        , L.length rFreqs
        ]
        1 .
      VU.fromList .
      -- L.map (\(a :+ b) -> (float2Double a) :+ (float2Double b)) .
      A.toList .
      runNWith
        ptx
        (gpuKernel''
           (A.constant ( sigma))
           (A.constant ( (log period)))
           freqArr) .
      A.fromList (A.Z A.:. (L.length xs)) -- .
      -- L.map
      --   (\(a, b, c, d, e) ->
      --      ( double2Float a
      --      , double2Float b
      --      , double2Float c
      --      , double2Float d
      --      , double2Float e)) 
        $
      xs

{-# INLINE computePinwheelCoefficients #-}
computePinwheelCoefficients ::
     Int
  -> Int
  -> Int
  -> Int
  -> PTX
  -> [(Double, Double, Double, Double, Complex Double)]
  -> Histogram (Complex Double)
computePinwheelCoefficients !maxPhiFreq !maxRhoFreq !maxThetaFreq !maxRFreq !ptx !xs =
  let freqFunc maxFreq = [-(fromIntegral maxFreq) .. (fromIntegral maxFreq)]
      phiFreqs = freqFunc maxPhiFreq
      rhoFreqs = freqFunc maxRhoFreq
      thetaFreqs = freqFunc maxThetaFreq
      rFreqs = freqFunc maxRFreq
      !freqArr =
        A.use $ computeFrequencyArray phiFreqs rhoFreqs thetaFreqs rFreqs
  in Histogram
       [ L.length phiFreqs
       , L.length rhoFreqs
       , L.length thetaFreqs
       , L.length rFreqs
       ]
       1 .
     VU.fromList .
     A.toList .
     runNWith ptx (pinwheelCoefficientsAcc freqArr) .
     A.fromList (A.Z A.:. (L.length xs)) $
     xs
