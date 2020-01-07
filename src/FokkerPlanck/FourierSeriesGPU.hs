{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ViewPatterns     #-}
module FokkerPlanck.FourierSeriesGPU where

import           Data.Array.Accelerate              as A
import qualified Data.Array.Accelerate.Data.Complex as AC
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Complex                       as C
import           Data.DList                         as DL
import           Data.List
import           Data.List                          as L
import           Data.Vector.Unboxed                as VU
import           FokkerPlanck.BrownianMotion        (Particle(..))
import           FokkerPlanck.FourierSeries
import           FokkerPlanck.Histogram
import           GHC.Float
import           Utils.Parallel

{-# INLINE rhoCutoff #-}
rhoCutoff :: (A.Ord a, A.Fractional a, A.Num a, Elt a) => Exp a -> Exp a
rhoCutoff !rho = cond (rho A.< 0) 0 rho

{-# INLINE moveParticle #-}  
moveParticle ::
     ( A.Ord a
     , A.Floating a
     , A.RealFloat a
     , A.Num a
     , Elt a
     , A.Fractional a
     )
  => Exp a
  -> Exp (a, a, a, a)
  -> Exp (a, a, a, a)
moveParticle !maxScale (unlift -> (phi, rho, theta, r)) =
  let !x = theta - phi
      !cosX = A.cos x
      !newPhi = phi + (A.atan2 (r * A.sin x) (rho + r * cosX))
      !newRho = A.sqrt $ (rho * rho + r * r + 2 * r * rho * cosX)
  in lift (newPhi, newRho, theta, r)

{-# INLINE coefficient #-}
coefficient ::
     (A.Floating a, A.Num a, A.RealFloat a, A.Elt (Complex a))
  => Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp (a, a, a, a)
  -> Exp (AC.Complex a)
coefficient halfLogPeriod rFreq thetaFreq rhoFreq phiFreq (unlift -> (phi, rho, theta, r)) =
  (lift (A.cos (phiFreq * phi + thetaFreq * (theta - phi)) AC.:+ 0)) *
  (AC.cis $
   (-A.pi) * (rhoFreq * (A.log rho) + rFreq * (A.log (r / rho))) /
   halfLogPeriod)

{-# INLINE computeFourierCoefficientsGPU #-}
computeFourierCoefficientsGPU ::
     Acc (A.Vector (Double, Double, Double, Double))
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> Exp Double
  -> Exp Double
  -> Exp (AC.Complex Double)
  -> PTX
  -> [DList Particle]
  -> Histogram (C.Complex Double)
computeFourierCoefficientsGPU !freqArr !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !deltaLogRho !maxScaleExp !halfLogPeriodExp !deltaLogRhoComplexExp !ptx !xs =
  let !particles = DL.toList . DL.concat $ xs
      !len = L.length particles
      !scaleSampledParticles =
        L.concat .
        parMap
          rdeepseq
          (\particle@(Particle phi rho theta r) ->
             let !n = Prelude.floor $ (log r + halfLogPeriod) / deltaLogRho
             in (phi, rho, theta, r) :
                [ ( phi
                  , rho
                  , theta
                  , (exp $
                     (Prelude.fromIntegral i * deltaLogRho - halfLogPeriod)))
                | i <- [0 .. n - 1]
                ]) $
        particles
      !movedParticleArr =
        compute .
        afst .
        A.filter
          (\particle ->
             let (_, rho, _, _) =
                   unlift particle :: ( Exp Double
                                      , Exp Double
                                      , Exp Double
                                      , Exp Double)
             in rho A./= 0) .
        A.map (moveParticle maxScaleExp) .
        use . A.fromList (Z :. (L.length scaleSampledParticles)) $
        scaleSampledParticles
  in Histogram
       [ L.length phiFreqs
       , L.length rhoFreqs
       , L.length thetaFreqs
       , L.length rFreqs
       ]
       len .
     VU.fromList .
     A.toList .
     runWith ptx .
     A.map
       (\(unlift -> (rFreq, thetaFreq, rhoFreq, phiFreq)) ->
          sfoldl
            (\s particle ->
               s +
               deltaLogRhoComplexExp *
               (coefficient
                  halfLogPeriodExp
                  rFreq
                  thetaFreq
                  rhoFreq
                  phiFreq
                  particle))
            0
            (constant Z)
            movedParticleArr) $
     freqArr
     


{-# INLINE computeFrequencyArray #-}
computeFrequencyArray ::
     [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Acc (A.Vector (Double, Double, Double, Double))
computeFrequencyArray !phiFreqs !rhoFreqs !thetaFreqs !rFreqs =
  use $
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
