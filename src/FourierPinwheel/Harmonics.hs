{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE StrictData       #-}
module FourierPinwheel.Harmonics where

import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List
import           Data.Vector.Generic as VG
import           Data.Vector.Unboxed as VU
import           Math.Gamma
import           Utils.Parallel

data FPData vector = FPData
  { getFPDataHarmonics       :: [vector]
  , getFPDataHarmonicsOffset :: [vector]
  , getFPDataCoef            :: vector
  , getFPDataCoefHollow      :: vector
  }

createHarmonics ::
     ( RealFloat e
     , Unbox e
     , NFData e
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     , Gamma (Complex e)
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> R.Array U DIM4 (Complex e)
  -> FPData (vector (Complex e))
createHarmonics numR2Freqs phiFreq rhoFreq thetaFreq rFreq alpha sigma periodR2 periodEnv coefficients =
  let radialConst = 2 * pi / (log periodEnv)
      centerR2 = div numR2Freqs 2
      harmonicVecs =
        parMap
          rdeepseq
          (\radialFreq ->
             createVector
               numR2Freqs
               (2 * phiFreq + 1)
               (pinwheel periodR2 radialConst alpha radialFreq))
          [-rhoFreq .. rhoFreq]
      harmonicVecsOffset =
        parMap
          rdeepseq
          (\radialFreq ->
             createVector
               numR2Freqs
               (2 * phiFreq + 1)
               (logpolarHarmonics periodR2 radialConst radialFreq))
          [rFreq,(rFreq - 1) .. -rFreq]
      coefficients2 = VU.convert . toUnboxed . computeS $ sumArr *^ coefficients
  in FPData harmonicVecs harmonicVecsOffset undefined coefficients2
  where
    sumArr =
      fromListUnboxed (extent coefficients) .
      parMap
        rdeepseq
        (\(r, theta, rho, phi) ->
           let angularFreq = phi - theta
               radialFreq = rho - r
               arr =
                 fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
                   let x = fromIntegral $ i - centerR2
                       y = fromIntegral $ j - centerR2
                       rho = sqrt $ x ^ 2 + y ^ 2
                       phi = atan2 y x
                   in if rho == 0
                        then 0
                        else pinwheel
                               periodR2
                               radialConst
                               alpha
                               radialFreq
                               angularFreq
                               phi
                               rho *
                             phaseShift angularFreq radialFreq radialConst alpha
           in (sumAllS arr) / (fromIntegral numR2Freqs ^ 2)) $
      [ (r, theta, rho, phi)
      | r <- [-rFreq .. rFreq]
      , theta <- [-thetaFreq .. thetaFreq]
      , rho <- [-rhoFreq .. rhoFreq]
      , phi <- [-phiFreq .. phiFreq]
      ]

{-# INLINE createVector #-}  
createVector ::
     (RealFloat e, Unbox e, VG.Vector vector (Complex e))
  => Int
  -> Int
  -> (Int -> e -> e -> Complex e)
  -> vector (Complex e)
createVector numR2Freqs numAngularFreqs f =
  let centerR2 = div numR2Freqs 2
      centerAngular = div numAngularFreqs 2
  in VG.convert .
     toUnboxed .
     computeS . fromFunction (Z :. numAngularFreqs :. numR2Freqs :. numR2Freqs) $ \(Z :. k :. i :. j) ->
       let x = fromIntegral $ i - centerR2
           y = fromIntegral $ j - centerR2
           rho = sqrt $ x ^ 2 + y ^ 2
           phi = atan2 y x
       in if rho == 0
            then 0
            else f (k - centerAngular) phi rho 

{-# INLINE pinwheel #-}
pinwheel :: (RealFloat e) => e -> e -> e -> Int -> Int -> e -> e -> Complex e
pinwheel periodR2 radialConst alpha radialFreq angularFreq phi rho =
  pi * (cis (fromIntegral angularFreq * phi)) *
  ((periodR2 / (pi * rho) :+ 0) **
   ((2 + alpha) :+ (radialConst * fromIntegral radialFreq)))

{-# INLINE logpolarHarmonics #-}
logpolarHarmonics ::
     (RealFloat e) => e -> e -> Int -> Int -> e -> e -> Complex e
logpolarHarmonics periodR2 radialConst radialFreq angularFreq phi rho =
  (cis (fromIntegral angularFreq * phi)) *
  ((periodR2 / (pi * rho) :+ 0) **
   (0 :+ (radialConst * fromIntegral radialFreq)))

{-# INLINE phaseShift #-}
phaseShift ::
     (RealFloat e, Gamma (Complex e)) => Int -> Int -> e -> e -> Complex e
phaseShift angularFreq radialFreq radialConst alpha =
  ((0 :+ (-1)) ^ (abs angularFreq)) *
  (gamma $
   ((2 + fromIntegral (abs angularFreq) + alpha) :+
    (radialConst * fromIntegral radialFreq)) /
   2) /
  (gamma $
   ((fromIntegral (abs angularFreq) - alpha) :+
    (radialConst * fromIntegral (-radialFreq))) /
   2)
