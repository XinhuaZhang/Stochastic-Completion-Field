{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE StrictData       #-}
module FourierPinwheel.Harmonics where

import           Data.Array.IArray   as IA
import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List
import           Data.List           as L
import           Data.Vector.Generic as VG
import           Data.Vector.Unboxed as VU
import           Filter.Utils
import           Math.Gamma
import           Utils.Parallel
import           Utils.List

data FPData vector = FPData
  { getFPDataHarmonics       :: [vector]
  , getFPDataHarmonicsOffset :: [vector]
  , getFPDataCoef            :: [vector]
  , getFPDataCoefHollow      :: [vector]
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
  -> R.Array U DIM4 (Complex e)
  -> IO (FPData (vector (Complex e)))
createHarmonics numR2Freqs phiFreq rhoFreq thetaFreq rFreq alpha periodR2 coefficients = do
  let harmonicVecs =
        parMap
          rdeepseq
          (\radialFreq ->
             createVector
               numR2Freqs
               phiFreq
               1
               (pinwheel periodR2 alpha radialFreq)) .
        getListFromNumber $
        rhoFreq
      harmonicVecsOffset =
        parMap
          rdeepseq
          (\radialFreq ->
             createVector
               numR2Freqs
               thetaFreq
               (-1)
               (logpolarHarmonics periodR2 radialFreq)) .
        L.map (\x -> -x) . getListFromNumber $
        rFreq
      centerPhiFreq = div phiFreq 2
      centerRhoFreq = div rhoFreq 2
      centerThetaFreq = div thetaFreq 2
      centerRFreq = div rFreq 2
  phaseShiftedCoefArr <-
    computeUnboxedP .
    R.zipWith (*) coefficients . fromFunction (extent coefficients) $ \(Z :. r :. theta :. rho :. phi) ->
      phaseShift
        (phi - centerPhiFreq - (theta - centerThetaFreq))
        (rho - centerRhoFreq - (r - centerRFreq))
        alpha
  let phaseShiftedCoefVecs =
        parMap
          rdeepseq
          (\r ->
             VU.convert . toUnboxed . computeS . R.slice phaseShiftedCoefArr $
             (Z :. r :. All :. All :. All))
          [0 .. rFreq - 1]
  coefHollowArr <- computeUnboxedP $ sumArr *^ coefficients
  let coefHollowVecs =
        parMap
          rdeepseq
          (\r ->
             VU.convert . toUnboxed . computeS . R.slice coefHollowArr $
             (Z :. r :. All :. All :. All))
          [0 .. rFreq - 1]
  return $
    FPData harmonicVecs harmonicVecsOffset phaseShiftedCoefVecs coefHollowVecs
  where
    !sumArr1 =
      computePinwheelArray
        numR2Freqs
        phiFreq
        rhoFreq
        thetaFreq
        rFreq
        alpha
        periodR2
    sumArr =
      fromListUnboxed (extent coefficients) .
      parMap
        rdeepseq
        (\(r, theta, rho, phi) ->
           let angularFreq = phi - theta
               radialFreq = rho - r
            in sumArr1 IA.! (radialFreq, angularFreq)) $
      [ (r, theta, rho, phi)
      | r <- getListFromNumber rFreq
      , theta <- getListFromNumber thetaFreq
      , rho <- getListFromNumber rhoFreq
      , phi <- getListFromNumber phiFreq
      ]

{-# INLINE createVector #-}
createVector ::
     (RealFloat e, Unbox e, VG.Vector vector (Complex e))
  => Int
  -> Int
  -> Int
  -> (Int -> e -> e -> Complex e)
  -> vector (Complex e)
createVector numR2Freqs numAngularFreqs sign f =
  let centerR2 = div numR2Freqs 2
      centerAngular = div numAngularFreqs 2
   in VG.convert .
      toUnboxed .
      computeS . 
      makeFilter2D . 
      fromFunction (Z :. numAngularFreqs :. numR2Freqs :. numR2Freqs) $ \(Z :. k :. i :. j) ->
        let x = fromIntegral (i - centerR2)
            y = fromIntegral (j - centerR2)
            rho = sqrt $ x ^ 2 + y ^ 2
            phi = atan2 y x
         in if rho == 0
              then 0
              else f (sign * (k - centerAngular)) phi rho

{-# INLINE pinwheel #-}
pinwheel :: (RealFloat e) =>  e -> e -> Int -> Int -> e -> e -> Complex e
pinwheel periodR2 alpha radialFreq angularFreq phi rho =
  pi * cis (fromIntegral angularFreq * phi) *
  ((periodR2 / (pi * rho) :+ 0) ** ((2 + alpha) :+ fromIntegral radialFreq))

{-# INLINE logpolarHarmonics #-}
logpolarHarmonics :: (RealFloat e) => e -> Int -> Int -> e -> e -> Complex e
logpolarHarmonics periodR2 radialFreq angularFreq phi rho =
  cis (fromIntegral angularFreq * phi) *
  ((periodR2 / (pi * rho) :+ 0) ** (0 :+ fromIntegral radialFreq))

{-# INLINE phaseShift #-}
phaseShift ::
     (RealFloat e, Gamma (Complex e)) => Int -> Int -> e -> Complex e
phaseShift angularFreq radialFreq alpha =
  ((0 :+ (-1)) ^ abs angularFreq) *
  gamma
    (((2 + fromIntegral (abs angularFreq) + alpha) :+ fromIntegral radialFreq) /
     2) /
  gamma
    (((fromIntegral (abs angularFreq) - alpha) :+ fromIntegral (-radialFreq)) /
     2)

{-# INLINE computeIndicesFromRadius #-}
computeIndicesFromRadius :: Double -> [(Int,Int)]
computeIndicesFromRadius radius =
  let rho = round radius
      r2 = rho ^ 2
   in L.filter
        (\(x, y) ->
           let r = fromIntegral $ x ^ 2 + y ^ 2
            in r < r2)
        [(x, y) | x <- [-rho .. rho], y <- [-rho .. rho]]

{-# INLINE computePinwheelArray #-}
computePinwheelArray ::
     ( Num e
     , RealFloat e
     , Gamma (Complex e)
     , Unbox e
     , NFData e
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> IA.Array (Int, Int) (Complex e)
computePinwheelArray numR2Freqs phiFreq rhoFreq thetaFreq rFreq alpha periodR2 =
  let (minPhiFreq, maxPhiFreq) = rangeFunc phiFreq
      (minRhoFreq, maxRhoFreq) = rangeFunc rhoFreq 
      (minThetaFreq, maxThetaFreq) = rangeFunc thetaFreq
      (minRFreq, maxRFreq) = rangeFunc rFreq      
      maxAngularFreq = max (maxPhiFreq + maxThetaFreq) (maxPhiFreq - minThetaFreq)
      minAngularFreq = min (minPhiFreq + minThetaFreq) (minPhiFreq - maxThetaFreq)
      maxRadialFreq = max (maxRhoFreq + maxRFreq) (maxRhoFreq - minRFreq)
      minRadialFreq = min (minRhoFreq + minRFreq) (minRhoFreq - maxRFreq)
      centerR2 = div numR2Freqs 2
      idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- [minRadialFreq .. maxRadialFreq]
        , angularFreq <- [minAngularFreq .. maxAngularFreq]
        ]
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             let arr =
                   fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
                     let x = fromIntegral $ i - centerR2
                         y = fromIntegral $ j - centerR2
                         rho = sqrt $ x ^ 2 + y ^ 2
                         phi = atan2 y x
                      in if rho == 0
                           then 0
                           else pinwheel
                                  periodR2
                                  alpha
                                  radialFreq
                                  angularFreq
                                  phi
                                  rho *
                                phaseShift angularFreq radialFreq alpha
              in sumAllS arr / (fromIntegral numR2Freqs ^ 2))
          idxs
   in listArray
        ((minRadialFreq, minAngularFreq), (maxRadialFreq, maxAngularFreq))
        pinwheels
  where
    rangeFunc n =
      let m = div n 2
       in if odd n
            then (-m, m)
            else (-m, m - 1)
