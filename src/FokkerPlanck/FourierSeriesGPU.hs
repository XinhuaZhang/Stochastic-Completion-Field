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

-- {-# INLINE computeFourierCoefficientsGPU #-}
computeFourierCoefficientsGPU ::
     GenIO
  -> Double
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> Acc (A.Vector (Double, Double, Double, Double))
  -> A.Exp Double
  -> A.Exp Double
  -> A.Exp (A.Complex Double)
  -> PTX
  -> [DList Particle]
  -> IO (Histogram (C.Complex Double))
computeFourierCoefficientsGPU !randomGen !thetaSigma !scaleSigma !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !deltaLogRho !freqArr !maxScaleExp !halfLogPeriodExp !deltaLogRhoComplexExp !ptx !xs = do
  let !particles = DL.toList . DL.concat $ xs
      -- zs = L.filter (\(_, rho, _, _) -> rho Prelude.> 0) . L.map (\(Particle a b c d) -> (a, b, c, d))  $ particles
  zs <-
    L.filter (\(_, rho, _, r) -> rho Prelude.> 0 && r Prelude.> 0) . L.concat <$>
    MP.mapM
       (\particle@(Particle phi rho theta r) -> do
          let !n = Prelude.floor $ r / deltaLogRho
              !delta = deltaLogRho / r
              !x =
                (\(Particle a b c d) -> (a, b, c, d)) .
                FokkerPlanck.BrownianMotion.moveParticle $
                particle
          gen <- createSystemRandom
          ys <-
            M.mapM
              (\i -> do
                 let thetaDist =
                       normalDistr 0 (thetaSigma * sqrt (delta * fromIntegral i))
                     rDist =
                       normalDistr 0 (scaleSigma * sqrt (delta * fromIntegral i))
                 deltaTheta <- genContVar thetaDist gen
                 deltaR <- genContVar rDist gen
                 return .
                   (\(Particle a b c d) ->
                      (a, b, c `thetaPlus` deltaTheta, r `scalePlus` deltaR)) .
                   FokkerPlanck.BrownianMotion.moveParticle $
                   (Particle phi rho theta (Prelude.fromIntegral i * deltaLogRho)))
              [1 .. n - 1]
          return (x : ys))
       particles
  let !len = L.length zs
  return .
    Histogram
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
    A.fromList (A.Z A.:. len) $
    zs


{-# INLINE computeFrequencyArray #-}
computeFrequencyArray ::
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> A.Vector (Double, Double, Double, Double)
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
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> PTX
  -> [(Double, Double, Double, Double, Double)]
  -> Histogram (Complex Double)
computeFourierCoefficientsGPU' !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !ptx !xs =
  let !freqArr =
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
     runNWith ptx (gpuKernel' freqArr) . A.fromList (A.Z A.:. (L.length xs)) $
     xs
