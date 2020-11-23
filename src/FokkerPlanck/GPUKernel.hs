{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ViewPatterns        #-}
module FokkerPlanck.GPUKernel where

import           Data.Array.Accelerate              as A
import           Data.Array.Accelerate.Data.Complex as A
import qualified Data.Complex                       as C
import Debug.Trace

{-# INLINE moveParticle #-}
moveParticle ::
     (A.RealFloat a, A.Num a, Elt a-- , FromIntegral Int a
     )
  => Exp (a, a, a, a)
  -> Exp (a, a, a, a)
moveParticle (unlift -> (phi, rho, theta, r)) =
  let !x = theta - phi
      !cosX = A.cos x
      !newPhi = phi + (A.atan2 (r * A.sin x) (rho + r * cosX))
      !newRho = A.sqrt $ (rho * rho + r * r + 2 * r * rho * cosX)
      -- !xx = A.round $ newRho * A.cos newPhi :: Exp Int
      -- !y = A.round $ newRho * A.sin newPhi :: Exp Int
      -- !newRho' = A.sqrt . A.fromIntegral $ xx * xx + y * y
      -- !newPhi' = A.atan2 (A.fromIntegral y) (A.fromIntegral xx)
  in lift (newPhi, newRho, theta, r)

{-# INLINE normalizationTerm #-}
normalizationTerm ::
     (A.Floating a, A.Num a, A.RealFloat a, A.Elt (Complex a))
  => Exp a
  -> Exp a
  -> Exp (Complex a)
normalizationTerm halfLogPeriod freq =
  (lift $ 1 A.:+ ((-A.pi) * freq / halfLogPeriod)) /
  (A.exp (lift $ halfLogPeriod A.:+ (-A.pi * freq)) -
   A.exp (lift $ (-halfLogPeriod) A.:+ (A.pi * freq)))

-- {-# INLINE coefficient #-}
-- coefficient ::
--      forall a. (A.Floating a, A.Num a, A.RealFloat a, A.Elt (Complex a), Prelude.Fractional a)
--   => Exp a
--   -> Exp a
--   -> Exp a
--   -> Exp a
--   -> Exp a
--   -> Exp (a, a, a, a)
--   -> Exp (A.Complex a)
-- coefficient halfLogPeriod rFreq thetaFreq rhoFreq phiFreq particle -- (unlift -> (phi, rho, theta, r))
--  =
--   let (phi, rho, theta, r) = unlift particle :: (Exp a, Exp a, Exp a, Exp a)
--      -- (normalizationTerm halfLogPeriod rhoFreq) *
--      -- (normalizationTerm halfLogPeriod rFreq) *
--      -- (lift (A.cos (phiFreq * phi + thetaFreq * (theta - phi)) A.:+ 0)) *
--                            -- (A.cis $
--                            --  (-2 * A.pi) * (rhoFreq * rho + rFreq * (r - rho)) /
--                            --  (A.exp halfLogPeriod))
     
--   in (lift $
--       ((A.exp $ (-0.5) * (rho + r)) *
--        (A.cos (phiFreq * phi  + thetaFreq * (theta - phi)))
--       ) A.:+
--       0) *
--      (A.cis $ (-1) * ((rhoFreq * rho + rFreq * (r - rho)
--                       ) -- * (2 * A.pi) / (A.exp halfLogPeriod) 
--                       -- + phiFreq * phi  + thetaFreq * (theta - phi)
--                      ))
--      -- (lift $ ((A.cos (phiFreq * phi + thetaFreq * (theta - phi) / 2))) A.:+ 0) *
--      -- (A.cis $ (-2 * A.pi) * (rhoFreq * rho + rFreq * (r - rho)) / (A.exp halfLogPeriod))
--      -- ((lift $ rho A.:+ 0) A.** (lift $ (-0.5) A.:+ (rFreq - rhoFreq))) *
--      -- ((lift $ r A.:+ 0) A.** (lift $ (-0.5) A.:+ (-rFreq)))
--      -- (A.cis $ (-A.pi) * (rhoFreq * rho + rFreq * (r - rho)
--      --                        ) / (halfLogPeriod))
--      -- (A.cis $
--      --  (-A.pi) * (rhoFreq * (A.log rho) + rFreq * (A.log (r / rho))) /
--      --  halfLogPeriod)
--   -- in ((lift $ rho A.:+ 0) A.**
--   --     (lift $ 0 A.:+ (-A.pi) * (rhoFreq - rFreq) / halfLogPeriod)) *
--   --    ((lift $ r A.:+ 0) A.** (lift $ 0 A.:+ (-A.pi) * rFreq / halfLogPeriod)) *
--   --    (A.cis $ (-phiFreq) * phi - thetaFreq * (theta - phi))
  
coefficient ::
     forall a. (A.Floating a, A.Num a, A.RealFloat a, A.Elt (Complex a), Prelude.Fractional a)
  => Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp (a, a, a, a)
  -> Exp (A.Complex a)
coefficient halfLogPeriod rFreq thetaFreq rhoFreq phiFreq particle =
  let (phi, rho, theta, r) = unlift particle :: (Exp a, Exp a, Exp a, Exp a)
  in (lift $
      ((A.exp $ (-0.5) * (rho + r)) *
       (A.cos (phiFreq * phi + thetaFreq * (theta - phi)))) A.:+
      0) *
     (A.cis $ (-1) * (rhoFreq * rho + rFreq * (r + rho)))

-- {-# INLINE gpuKernel #-}
gpuKernel ::
     forall a. (A.Eq a, A.Floating a, A.Num a, A.RealFloat a, A.Elt (Complex a), A.FromIntegral Int a, Prelude.Fractional a)
  => Exp a
  -> Exp a
  -> Exp (A.Complex a)
  -> Acc (A.Vector (a, a, a, a))
  -> Acc (A.Vector (a, a, a, a))
  -> Acc (A.Vector (A.Complex a))
gpuKernel !maxScaleExp !halfLogPeriodExp !deltaLogRhoComplexExp freqArr particles =
  let -- delta = constant 0.01
      movedParticleArr =
        compute .
        A.map
          (\particle ->
             let (phi, rho, theta, r) =
                   unlift particle :: (Exp a, Exp a, Exp a, Exp a)
                 logRho = A.log rho 
                   -- (A.fromIntegral $
                   --  (A.round (((A.log rho) + halfLogPeriodExp) / delta) :: Exp Int)) *
                   -- delta -
                   -- halfLogPeriodExp
                 -- newRho =
                 --   (A.fromIntegral $ (A.round (rho / delta) :: Exp Int)) * delta :: Exp a
             in lift (phi, logRho, theta, (A.log r) :: Exp a)) -- .
        -- afst .
        -- A.filter
        --   (\particle ->
        --      let (_, rho, _, _) =
        --            unlift particle :: (Exp a, Exp a, Exp a, Exp a)
        --      in rho A.> (A.constant 0) ) -- .
        -- A.map moveParticle 
        $
        particles
  in A.map
       (\(unlift -> (rFreq, thetaFreq, rhoFreq, phiFreq))
          -- (lift $ delta A.:+ 0) /
          -- (lift $ (8 * A.pi * A.pi * halfLogPeriodExp) A.:+ 0) *
         ->
          (sfoldl
             (\s particle ->
                s +
                (coefficient
                   halfLogPeriodExp
                   rFreq
                   thetaFreq
                   rhoFreq
                   phiFreq
                   particle))
             0
             (constant Z)
             movedParticleArr) * deltaLogRhoComplexExp
          -- (lift $
          --  (16 * A.pi * A.pi * halfLogPeriodExp * halfLogPeriodExp) A.:+ 0)
             -- movedParticleArr
        ) $
     freqArr

{-# INLINE convolveKernel #-}
convolveKernel ::
     (A.Num a, A.RealFloat a, A.Elt (Complex a))
  => Acc (Array DIM4 (Complex a))
  -> Acc (Array DIM4 (Complex a))
  -> Acc (Scalar Int)
  -> Acc (Scalar Int)
  -> Acc (Array DIM4 (Complex a))
  -> Acc (Vector (Complex a))
convolveKernel coefficients harmonics thetaIdx rIdx input =
  let (Z :. numRhoFreq :. numPhiFreq :. cols :. rows) =
        unlift . shape $ input :: (Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int)
      coefficientsArr =
        A.replicate (lift (Z :. All :. All :. cols :. rows)) .
        slice coefficients $
        (lift (Z :. (the rIdx) :. (the thetaIdx) :. All :. All))
      harmonicsArr =
        backpermute
          (shape input)
          (\(unlift -> Z :. rho :. phi :. col :. row :: (Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int)) ->
             lift (Z :. (rho + the rIdx) :. (phi + the thetaIdx) :. col :. row))
          harmonics
  in flatten .
     A.sum .
     A.sum .
     backpermute
       (lift (Z :. cols :. rows :. numRhoFreq :. numPhiFreq))
       (\(unlift -> Z :. col :. row :. rho :. phi :: (Z :. Exp Int :. Exp Int :. Exp Int :. Exp Int)) ->
          lift (Z :. rho :. phi :. col :. row)) .
     A.zipWith (*) harmonicsArr . A.zipWith (*) coefficientsArr $
     input


{-# INLINE coefficient' #-}
coefficient' ::
     forall a.
     ( A.Floating a
     , A.Num a
     , A.RealFloat a
     , A.Elt (Complex a)
     , Prelude.Fractional a
     )
  => Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp (a, a, a, a, a)
  -> Exp (A.Complex a)
coefficient' sigma rFreq thetaFreq rhoFreq phiFreq particle =
  let (phi, rho, theta, r, v) =
        unlift particle :: (Exp a, Exp a, Exp a, Exp a, Exp a)
  in (lift $ (v * (A.exp $ (sigma - 1) * (rho + r))) A.:+ 0) *
     (A.cis $
      (-1) *
      (rhoFreq * rho + rFreq * (r - rho) + phiFreq * phi +
       thetaFreq * (theta - phi))) 

gpuKernel' ::
     forall a.
     ( A.Eq a
     , A.Floating a
     , A.Num a
     , A.RealFloat a
     , A.Elt (Complex a)
     , A.FromIntegral Int a
     , Prelude.Fractional a
     )
  => Exp a
  -> Acc (A.Vector (a, a, a, a))
  -> Acc (A.Vector (a, a, a, a, a))
  -> Acc (A.Vector (A.Complex a))
gpuKernel' sigma freqArr xs =
  A.map
    (\(unlift -> (rFreq, thetaFreq, rhoFreq, phiFreq)) ->
       A.sfoldl
         (\s particle ->
            s + (coefficient' sigma rFreq thetaFreq rhoFreq phiFreq particle))
         0
         (constant Z)
         xs)
    freqArr
    

{-# INLINE coefficient'' #-}
coefficient'' ::
     forall a.
     ( A.Floating a
     , A.Num a
     , A.RealFloat a
     , A.Elt (Complex a)
     , Prelude.Fractional a
     )
  => Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp (a, a, a, a, a)
  -> Exp (A.Complex a)
coefficient'' sigma  rFreq thetaFreq rhoFreq phiFreq particle =
  let (phi, rho, theta, r, v) =
        unlift particle :: (Exp a, Exp a, Exp a, Exp a, Exp a)
   in lift
        ((v * A.exp ((sigma - 1) * rho) *
          A.cos (phiFreq * phi + thetaFreq * (theta - phi))
         ) :+
         0) *
      A.cis ((-1) * ((rhoFreq * rho + rFreq * (r - rho)))) 


gpuKernel'' ::
     forall a.
     ( A.Eq a
     , A.Floating a
     , A.Num a
     , A.RealFloat a
     , A.Elt (Complex a)
     , A.FromIntegral Int a
     , Prelude.Fractional a
     )
  => Exp a
  -> Acc (A.Vector (a, a, a, a))
  -> Acc (A.Vector (a, a, a, a, a))
  -> Acc (A.Vector (A.Complex a))
gpuKernel'' sigma freqArr xs =
  A.map
    (\(unlift -> (rFreq, thetaFreq, rhoFreq, phiFreq)) ->
       A.sfoldl
         (\s particle ->
            s +
            coefficient'' sigma rFreq thetaFreq rhoFreq phiFreq particle)
         0
         (constant Z)
         xs)
    freqArr
    

{-# INLINE pinwheelAcc #-}
pinwheelAcc ::
     forall a. (A.Floating a, A.Num a, A.RealFloat a, A.Elt (Complex a), Prelude.Fractional a)
  => Exp a
  -> Exp a
  -> Exp a
  -> Exp a
  -> Exp (a, a, a, a, A.Complex a)
  -> Exp (A.Complex a)
pinwheelAcc rFreq thetaFreq rhoFreq phiFreq particle =
  let (phi, rho, theta, r, v) =
        unlift particle :: (Exp a, Exp a, Exp a, Exp a, Exp (A.Complex a))
  in v * (lift $ ((A.exp $ (-0.5) * (rho + r))) A.:+ 0) *
     (A.cis $
      (-1) *
      (rhoFreq * rho + rFreq * (r - rho) + phiFreq * phi +
       thetaFreq * (theta - phi))) 

pinwheelCoefficientsAcc ::
     forall a.
     ( A.Eq a
     , A.Floating a
     , A.Num a
     , A.RealFloat a
     , A.Elt (Complex a)
     , A.FromIntegral Int a
     , Prelude.Fractional a
     )
  => Acc (A.Vector (a, a, a, a))
  -> Acc (A.Vector (a, a, a, a, A.Complex a))
  -> Acc (A.Vector (A.Complex a))
pinwheelCoefficientsAcc freqArr xs =
  A.map
    (\(unlift -> (rFreq, thetaFreq, rhoFreq, phiFreq)) ->
       A.sfoldl
         (\s particle ->
            s + (pinwheelAcc rFreq thetaFreq rhoFreq phiFreq particle))
         0
         (constant Z)
         xs)
    freqArr
