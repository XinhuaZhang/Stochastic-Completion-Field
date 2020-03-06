{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
module FokkerPlanck.FourierSeriesGPU where

import           Data.Array.Accelerate              as A
import qualified Data.Array.Accelerate.Data.Complex as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Complex                       as C
import           Data.DList                         as DL
import           Data.List                          as L
import           Data.Vector.Unboxed                as VU
import           FokkerPlanck.BrownianMotion        (Particle (..),moveParticle)
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.GPUKernel
import           FokkerPlanck.Histogram
import           GHC.Float
import           Utils.Parallel

{-# INLINE computeFourierCoefficientsGPU #-}
computeFourierCoefficientsGPU ::
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> Acc (A.Vector (Double, Double, Double, Double))
  -> Exp Double
  -> Exp Double
  -> Exp (A.Complex Double)
  -> PTX
  -> [DList Particle]
  -> Histogram (C.Complex Double)
computeFourierCoefficientsGPU !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !deltaLogRho !freqArr !maxScaleExp !halfLogPeriodExp !deltaLogRhoComplexExp !ptx !xs =
  let !particles = DL.toList . DL.concat $ xs
      !len = L.length particles
      !scaleSampledParticles =
        L.concat .
        parMap
          rdeepseq
          (\particle@(Particle phi rho theta r) ->
             let n = Prelude.floor $ r / deltaLogRho
             in L.map (\(Particle a b c d) -> (a, b, c, r)) -- .
                -- L.map FokkerPlanck.BrownianMotion.moveParticle
             $
                (Particle phi rho theta r) :
                [ -- (Particle phi rho theta (Prelude.fromIntegral i * deltaLogRho))
                -- | i <- [1 .. n - 1]
                ]) $
        particles
      -- !scaleSampledParticles =
      --   L.concat .
      --   parMap
      --     rdeepseq
      --     (\particle@(Particle phi rho theta r) ->
      --        let !n = Prelude.floor $ (log r + halfLogPeriod) / deltaLogRho
      --        in (phi, rho, theta, r) :
      --           [ ( phi
      --             , rho
      --             , theta
      --             , (exp $
      --                (Prelude.fromIntegral i * deltaLogRho - halfLogPeriod)))
      --           | i <- [0 .. n - 1]
      --           ]) $
      --   particles
  in Histogram
       [ L.length phiFreqs
       , L.length rhoFreqs
       , L.length thetaFreqs
       , L.length rFreqs
       ]
       len .
     VU.fromList .
     A.toList .
     runNWith
       ptx
       (gpuKernel maxScaleExp halfLogPeriodExp deltaLogRhoComplexExp freqArr) .
     A.fromList (Z :. (L.length scaleSampledParticles)) $
     scaleSampledParticles

{-# INLINE computeFrequencyArray #-}
computeFrequencyArray ::
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> A.Vector (Double, Double, Double, Double)
computeFrequencyArray !phiFreqs !rhoFreqs !thetaFreqs !rFreqs =
  A.fromList
    (Z :.
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
