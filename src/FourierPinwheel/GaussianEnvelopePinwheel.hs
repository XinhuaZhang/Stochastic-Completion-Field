{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
module FourierPinwheel.GaussianEnvelopePinwheel where

import           Data.Array.Repa             as R
import           Data.Complex
import           Data.Vector.Generic         as VG
import           Data.Vector.Unboxed         as VU
import           FourierPinwheel.Hypergeo1F1
import           Math.Gamma
import           Pinwheel.FourierSeries2D
import           Utils.Distribution
import           Utils.Parallel

-- The Fourier coefficients of e^{-a^2 r^2} r^\mu
{-# INLINE gaussianPinwheelFourierCoefficients #-}
gaussianPinwheelFourierCoefficients ::
     (RealFloat a, Gamma (Complex a), Enum a)
  => Int
  -> a
  -> a
  -> a
  -> Int
  -> Int
  -> a
  -> a
  -> a
  -> Complex a
gaussianPinwheelFourierCoefficients numR2Freqs periodR2 a sigma angularFreq radialFreq periodEnv phi rho =
  let periodConst = (-2) * pi / log periodEnv
      piRhoPConst = pi * rho / periodR2
      real = pi / periodR2 ^ 2 * piRhoPConst ^ abs angularFreq
      mu =
        (2 + sigma + fromIntegral (abs angularFreq)) :+
        (periodConst * fromIntegral radialFreq)
      alpha = mu / 2
      beta = fromIntegral (1 + abs angularFreq) :+ 0
      z = ((-1) * (piRhoPConst / a) ^ 2) :+ 0
      img =
        (0 :+ (-1)) ^ abs angularFreq * cis (fromIntegral (-angularFreq) * phi) *
        gamma alpha /
        gamma beta *
        (a :+ 0) ** (-mu) *
        hypergeom alpha beta z
   in (real :+ 0) * img

gaussianPinwheel ::
     ( RealFloat a
     , Gamma (Complex a)
     , Enum a
     , Unbox a
     , VG.Vector vector (Complex a)
     , NFData (vector (Complex a))
     )
  => Int
  -> a
  -> a
  -> a
  -> Int
  -> Int
  -> a
  -> a
  -> a
  -> vector (Complex a)
gaussianPinwheel numR2Freqs periodR2 stdR2 sigma thetaFreq rFreq periodEnv stdTheta stdR =
  let zeroVec = VG.replicate (numR2Freqs ^ 2) 0
      a = 1 / (stdR2 * sqrt 2)
   in VG.concat .
      parMap
        rdeepseq
        (\(radialFreq, angularFreq) ->
           if angularFreq == 0
             then let pinwheel =
                        centerHollowArray numR2Freqs $
                        createFrequencyArray
                          numR2Freqs
                          (gaussianPinwheelFourierCoefficients
                             numR2Freqs
                             periodR2
                             a
                             sigma
                             angularFreq
                             radialFreq
                             periodEnv)
                      arr =
                        R.map
                          (* ((gaussian1DFourierCoefficients
                                 (fromIntegral radialFreq)
                                 (log periodEnv)
                                 stdR) :+
                              0)) $
                        centerHollowArray numR2Freqs pinwheel
                   in VG.convert . toUnboxed . computeS $ arr
             else zeroVec) $
      [ (radialFreq, angularFreq)
      | radialFreq <- [-rFreq .. rFreq]
      , angularFreq <- [-thetaFreq .. thetaFreq]
      ]
      

-- This one has a orientation preference at 0 degree. It is used for Koffka cross problem.
gaussianPinwheel1 ::
     ( RealFloat a
     , Gamma (Complex a)
     , Enum a
     , Unbox a
     , VG.Vector vector (Complex a)
     , NFData (vector (Complex a))
     )
  => Int
  -> a
  -> a
  -> a
  -> Int
  -> Int
  -> a
  -> a
  -> a
  -> vector (Complex a)
gaussianPinwheel1 numR2Freqs periodR2 stdR2 sigma thetaFreq rFreq periodEnv stdTheta stdR =
  let a = 1 / (stdR2 * sqrt 2)
   in VG.concat .
      parMap
        rdeepseq
        (\(radialFreq, angularFreq) ->
           let pinwheel =
                 centerHollowArray numR2Freqs $
                 createFrequencyArray
                   numR2Freqs
                   (gaussianPinwheelFourierCoefficients
                      numR2Freqs
                      periodR2
                      a
                      sigma
                      0
                      radialFreq
                      periodEnv)
               arr =
                 R.map
                   (* ((gaussian1DFreq (fromIntegral angularFreq) stdTheta *
                        gaussian1DFourierCoefficients
                          (fromIntegral radialFreq)
                          (log periodEnv)
                          stdR) :+
                       0)) $
                 centerHollowArray numR2Freqs pinwheel
            in VG.convert . toUnboxed . computeS $ arr) $
      [ (radialFreq, angularFreq)
      | radialFreq <- [-rFreq .. rFreq]
      , angularFreq <- [-thetaFreq .. thetaFreq]
      ]
