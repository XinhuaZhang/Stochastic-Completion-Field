{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE Strict              #-}
module Pinwheel.FourierSeries2D where

import           Control.Monad                  as M
import           Control.Monad.IO.Class
import           Control.Monad.Trans.Resource
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX as A
import           Data.Array.IArray              as IA
import           Data.Array.Repa                as R
import           Data.Complex
import           Data.Conduit                   as C
import           Data.Conduit.List              as CL
import           Data.List                      as L
import           Data.Vector.Generic            as VG
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.FourierSeries2D
import           Math.Gamma
import           Pinwheel.Base
import           Pinwheel.List
import           Utils.BLAS
import           Utils.Distribution
import           Utils.List
import           Utils.Parallel                 hiding ((.|))
import           Utils.SimpsonRule
import           Utils.Time
import DFT.Plan
import Filter.Utils
import Control.Concurrent.Async
import Numeric.GSL.Special.Bessel

pinwheelFourierSeries ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> Int
  -> Int
  -> IO (IA.Array (Int, Int) (vector (Complex e)))
pinwheelFourierSeries deviceIDs ptxs numR2Freqs numPoints delta periodR2 phiFreq rhoFreq thetaFreq rFreq sigma numBatchR2 numBatchPinwheelFreqs = do
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- pinwheelFreqs rhoFreq rFreq
        , angularFreq <- pinwheelFreqs phiFreq thetaFreq
        ]
      centerR2Freq = div numR2Freqs 2
  pinwheels <-
    computeUnboxedP .
    R.traverse
      (fromListUnboxed (Z :. (L.length idxs)) idxs)
      (\(Z :. freqs) -> (Z :. freqs :. numR2Freqs :. numR2Freqs)) $ \f (Z :. pFreq :. xFreq :. yFreq) ->
      fourierMellin
        sigma
        (snd $ f (Z :. pFreq))
        (fst $ f (Z :. pFreq))
        ( delta * fromIntegral (xFreq - centerR2Freq)
        , delta * fromIntegral (yFreq - centerR2Freq))
  -- pinwheels <-
  --   M.mapM
  --     (\idx ->
  --        fmap
  --          (CuMat (numR2Freqs ^ 2) (L.length idx) .
  --           CuVecHost . VU.convert . toUnboxed) .
  --        computeP .
  --        R.traverse
  --          (fromListUnboxed (Z :. (L.length idx)) idx)
  --          (\(Z :. freqs) -> (Z :. numR2Freqs :. numR2Freqs :. freqs)) $ \f (Z :. xFreq :. yFreq :. pFreq) ->
  --          fourierMellin
  --            sigma
  --            (snd $ f (Z :. pFreq))
  --            (fst $ f (Z :. pFreq))
  --            ( fromIntegral (xFreq - centerR2Freq)
  --            , fromIntegral (yFreq - centerR2Freq))) .
  --   divideListN numBatchPinwheelFreqs $
  --   idxs
  -- printCurrentTime "Compute Fourier Series... "
  -- arr <-
  --   computeFourierSeriesR2Stream
  --     deviceIDs
  --     ptxs
  --     numR2Freqs
  --     numPoints
  --     periodR2
  --     delta
  --     numBatchR2
  --     pinwheels
  -- printCurrentTime "Compute Fourier Series Done"
  let (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound phiFreq thetaFreq
  return .
    listArray ((radialLB, angularLB), (radialUB, angularUB)) .
    parMap
      rdeepseq
      (\i ->
         VG.convert . toUnboxed . computeS . R.slice pinwheels $  -- arr $
         (Z :. i :. All :. All)) $
    [0 .. (L.length idxs - 1)]


-- Stream
{-# INLINE source #-}
source :: Int -> Int -> Int -> Int -> ConduitT () (Int, Int) (ResourceT IO) ()
source phiFreq rhoFreq thetaFreq rFreq =
  CL.sourceList
    [ (radialFreq, angularFreq)
    | radialFreq <- pinwheelFreqs rhoFreq rFreq
    , angularFreq <- pinwheelFreqs phiFreq thetaFreq
    ]

conduit ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , VG.Vector vector (Complex e)
     )
  => [Int]
  -> Int
  -> Int
  -> e
  -> Int
  -> e
  -> [[CuMat (Complex e)]]
  -> ConduitT (Int, Int) [vector (Complex e)] (ResourceT IO) ()
conduit deviceIDs numR2Freqs numPoints periodR2 batchSize sigma inverseR2Harmonics = do
  liftIO $ printCurrentTime "conduit"
  idx <- CL.take batchSize
  unless
    (L.null idx)
    (do let centerR2Freq = div numR2Freqs 2
        pinwheels <-
          fmap
            (CuMat (L.length idx) (numR2Freqs ^ 2) .
             CuVecHost . VU.convert . toUnboxed) .
          liftIO .
          computeP .
          R.traverse
            (fromListUnboxed (Z :. (L.length idx)) idx)
            (\(Z :. freqs) -> (Z :. freqs :. numR2Freqs :. numR2Freqs)) $ \f (Z :. pFreq :. xFreq :. yFreq) ->
            fourierMellin
              sigma
              (snd $ f (Z :. pFreq))
              (fst $ f (Z :. pFreq))
              ( fromIntegral (xFreq - centerR2Freq)
              , fromIntegral (yFreq - centerR2Freq))
        arr <-
          liftIO $
          computeFourierSeriesR2
            deviceIDs
            numR2Freqs
            numPoints
            periodR2
            inverseR2Harmonics
            [pinwheels]
        yield .
          L.map
            (\i ->
               VG.convert . toUnboxed . computeS . R.slice arr $
               (Z :. i :. All :. All)) $
          [0 .. (L.length idx - 1)])

pinwheelFourierSeriesStream ::
     ( Storable e
     , CUBLAS (Complex e)
     , Unbox e
     , RealFloat e
     , A.Elt e
     , A.Elt (Complex e)
     , Floating (A.Exp e)
     , A.FromIntegral Int e
     , VG.Vector vector (Complex e)
     )
  => [Int]
  -> [PTX]
  -> Int
  -> Int
  -> e
  -> e
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> Int
  -> Int
  -> Int
  -> IO (IA.Array (Int, Int) (vector (Complex e)))
pinwheelFourierSeriesStream deviceIDs ptxs numR2Freqs numPoints delta periodR2 phiFreq rhoFreq thetaFreq rFreq sigma numBatchR2 batchSize numBatchPinwheelFreqs = do
  inverseR2Harmonics <-
    createInverseHarmonicMatriesGPU
      ptxs
      numBatchR2
      numPoints
      numR2Freqs
      periodR2
      delta
  let (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound rhoFreq rFreq
  fmap (listArray ((radialLB, angularLB), (radialUB, radialUB)) . L.concat) .
    runConduitRes $
    source phiFreq rhoFreq thetaFreq rFreq .|
    conduit
      deviceIDs
      numR2Freqs
      numPoints
      periodR2
      batchSize
      sigma
      inverseR2Harmonics .|
    CL.consume

-- 2D spatial domain to frequency domain:
-- transform: cis (-freq * x)   
-- inverse transform: cis (freq * x)   

-- Analytical solution type 1
-- The pinwheel harmonics are cis (-freq * x)   
{-# INLINE analyticalFourierSeriesFunc1 #-}
analyticalFourierSeriesFunc1 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e))
  => Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> e
  -> Complex e
analyticalFourierSeriesFunc1 angularFreq radialFreq sigma periodR2 periodEnv phi rho =
  let radialConst = 2 * pi / log periodEnv
   in pi * ((0 :+ 1) ^ abs angularFreq) *
      cis (fromIntegral (-angularFreq) * phi) *
      ((periodR2 / (pi * rho) :+ 0) **
       ((2 + sigma) :+ (radialConst * fromIntegral (-radialFreq)))) *
      gamma
        (((2 + fromIntegral (abs angularFreq) + sigma) :+
          (radialConst * fromIntegral (-radialFreq))) /
         2) /
      gamma
        (((fromIntegral (abs angularFreq) - sigma) :+
          (radialConst * fromIntegral radialFreq)) /
         2) 

-- The pinwheel is in the frequency domain, the function computes the Fourier series
{-# INLINE analyticalFourierSeries1 #-}
analyticalFourierSeries1 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierSeries1 numPoints delta angularFreq radialFreq sigma periodR2 periodEnv =
  let center = div numPoints 2
  in fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
       let x = fromIntegral $ i - center
           y = fromIntegral $ j - center
           rho = sqrt $ x ^ 2 + y ^ 2
           phi = atan2 y x
       in if rho <= 0
            then 0
            else analyticalFourierSeriesFunc1
                   angularFreq
                   radialFreq
                   sigma
                   periodR2
                   periodEnv
                   phi
                   (rho * delta)

-- The pinwheel is in the spatial domain, the function computes the Fourier coefficients
{-# INLINE analyticalFourierCoefficients1 #-}
analyticalFourierCoefficients1 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierCoefficients1 numFreqs delta angularFreq radialFreq sigma periodR2 periodEnv =
  let c = ((-1) ^ (abs angularFreq)) :+ 0
  in R.map (* c) $
     analyticalFourierSeries1 numFreqs delta angularFreq radialFreq sigma periodR2 periodEnv


-- Analytical solution type 2
-- The pinwheel harmonics are cis (freq * x)   
{-# INLINE analyticalFourierSeriesFunc2 #-}
analyticalFourierSeriesFunc2 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e))
  => Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> e
  -> Complex e
analyticalFourierSeriesFunc2 angularFreq radialFreq sigma periodR2 periodEnv phi rho =
  let radialConst = 2 * pi / (log periodEnv)
  in pi * ((0 :+ 1) ^ (abs angularFreq)) *
     ((cis (fromIntegral angularFreq * phi))) *
     ((periodR2 / (pi * rho) :+ 0) **
      ((2 + sigma) :+ (radialConst * fromIntegral radialFreq))) *
     (gamma $
      ((2 + fromIntegral (abs angularFreq) + sigma) :+
       (radialConst * fromIntegral radialFreq)) /
      2) /
     (gamma $
      ((fromIntegral (abs angularFreq) - sigma) :+
       (radialConst * fromIntegral (-radialFreq))) /
      2) 
  

-- The pinwheel is in the frequency domain, the function computes the Fourier series
{-# INLINE analyticalFourierSeries2 #-}
analyticalFourierSeries2 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierSeries2 numPoints delta angularFreq radialFreq sigma periodR2 periodEnv =
  let center = div numPoints 2
      arr =
        fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
          let x = delta * (fromIntegral $ i - center)
              y = delta * (fromIntegral $ j - center)
              rho = sqrt $ x ^ 2 + y ^ 2
              phi = atan2 y x
          in if rho == 0 -- || (rho > period / 2) -- || phi /= 0
               then 0
               else analyticalFourierSeriesFunc2
                      angularFreq
                      radialFreq
                      sigma
                      periodR2
                      periodEnv
                      phi
                      rho
     -- if angularFreq == 0
     --   then fromFunction (Z :. numPoints :. numPoints) $ \_ -> 0
     --   else
  in arr

-- The pinwheel is in the spatial domain, the function computes the Fourier coefficients
{-# INLINE analyticalFourierCoefficients2 #-}
analyticalFourierCoefficients2 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierCoefficients2 numFreqs delta angularFreq radialFreq sigma periodR2 periodEnv =
  let c = ((-1) ^ (abs angularFreq)) / periodR2 :+ 0
  in R.map (* c) $
     analyticalFourierSeries2
       numFreqs
       delta
       angularFreq
       radialFreq
       sigma
       periodR2
       periodEnv
       

-- Analytical solution type 3
-- The pinwheel harmonics are cis (-freq * x)   
-- The envelope is r^2 * r^\alpha, where \alpha \in (-2,-0.5)
{-# INLINE analyticalFourierSeriesFunc3 #-}
analyticalFourierSeriesFunc3 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e))
  => Int
  -> Int
  -> e
  -> e
  -> e
  -> e
  -> e
  -> Complex e
analyticalFourierSeriesFunc3 angularFreq radialFreq sigma periodR2 periodEnv phi rho =
  let radialConst = 2 * pi / log periodEnv
   in ((fromIntegral angularFreq ^ 2 :+ 0) -
       ((sigma + 2) :+ (radialConst * fromIntegral (-radialFreq))) ^ 2) /
      (8 * pi) *
      ((0 :+ 1) ^ abs angularFreq) *
      cis (fromIntegral (-angularFreq) * phi) /
      (rho ^ 2 :+ 0) *
      ((periodR2 / (pi * rho) :+ 0) **
       ((2 + sigma) :+ (radialConst * fromIntegral (-radialFreq)))) *
      gamma
        (((2 + fromIntegral (abs angularFreq) + sigma) :+
          (radialConst * fromIntegral (-radialFreq))) /
         2) /
      gamma
        (((fromIntegral (abs angularFreq) - sigma) :+
          (radialConst * fromIntegral radialFreq)) /
         2) 
         
{-# INLINE analyticalFourierSeries3 #-}
analyticalFourierSeries3 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierSeries3 numPoints delta angularFreq radialFreq sigma periodR2 periodEnv =
  let center = div numPoints 2
  in fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
       let x = fromIntegral $ i - center
           y = fromIntegral $ j - center
           rho = sqrt $ x ^ 2 + y ^ 2
           phi = atan2 y x
       in if rho == 0
            then 0
            else analyticalFourierSeriesFunc3
                   angularFreq
                   radialFreq
                   sigma
                   periodR2
                   periodEnv
                   phi
                   (rho * delta)
         
{-# INLINE analyticalFourierCoefficients3 #-}
analyticalFourierCoefficients3 ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierCoefficients3 numFreqs delta angularFreq radialFreq sigma periodR2 periodEnv =
  let c = ((-1) ^ (abs angularFreq)) :+ 0
  in R.map (* c) $
     analyticalFourierSeries3 numFreqs delta angularFreq radialFreq sigma periodR2 periodEnv


pinwheelFourierCoefficientsAnatical ::
     ( Unbox e
     , RealFloat e
     , Gamma (Complex e)
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> IA.Array (Int, Int) (vector (Complex e))
pinwheelFourierCoefficientsAnatical numR2Freqs phiFreq rhoFreq thetaFreq rFreq sigma periodR2 periodEnv =
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- pinwheelFreqs rhoFreq rFreq
        , angularFreq <- pinwheelFreqs phiFreq thetaFreq
        ]
      centerR2Freq = div numR2Freqs 2
      (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound phiFreq thetaFreq
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             VG.convert .
             toUnboxed .
             computeS -- .
             -- R.zipWith
             --   (*)
             --   (fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
             --      gaussian2D
             --        (fromIntegral $ i - centerR2Freq)
             --        (fromIntegral $ j - centerR2Freq)
             --        (fromIntegral $ div centerR2Freq 1))
            $
             analyticalFourierCoefficients2
               numR2Freqs
               1
               (angularFreq)
               (radialFreq)
               sigma
               periodR2
               periodEnv)
          idxs
  in listArray ((radialLB, angularLB), (radialUB, angularUB)) pinwheels
  

pinwheelFourierCoefficientsAnaticalList ::
     ( Unbox e
     , RealFloat e
     , Gamma (Complex e)
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> [vector (Complex e)]
pinwheelFourierCoefficientsAnaticalList numR2Freqs phiFreq rhoFreq thetaFreq rFreq sigma periodR2 periodEnv =
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- [-rhoFreq .. rhoFreq]
        , angularFreq <- [-phiFreq .. phiFreq]
        ]
      centerR2Freq = div numR2Freqs 2
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             VG.convert . toUnboxed . computeS $
             analyticalFourierCoefficients2
               numR2Freqs
               1
               (angularFreq)
               (radialFreq)
               sigma
               periodR2
               periodEnv)
          idxs
  in pinwheels
  
pinwheelFourierCoefficientsAnaticalList1 ::
     ( Unbox e
     , RealFloat e
     , Gamma (Complex e)
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> [vector (Complex e)]
pinwheelFourierCoefficientsAnaticalList1 numR2Freqs phiFreq rhoFreq thetaFreq rFreq sigma periodR2 periodEnv =
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- [-rhoFreq .. rhoFreq]
        , angularFreq <- [-phiFreq .. phiFreq]
        ]
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             VG.convert . toUnboxed . computeS $
             analyticalFourierCoefficients2
               numR2Freqs
               1
               (angularFreq)
               (radialFreq)
               sigma
               periodR2
               periodEnv)
          idxs
  in pinwheels
  
pinwheelFourierCoefficientsAnaticalList2 ::
     ( Unbox e
     , RealFloat e
     , Gamma (Complex e)
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> [vector (Complex e)]
pinwheelFourierCoefficientsAnaticalList2 numR2Freqs phiFreq rhoFreq thetaFreq rFreq sigma periodR2 periodEnv =
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- [rFreq,(rFreq - 1) .. -rFreq]
        , angularFreq <- [thetaFreq,(thetaFreq - 1) .. -thetaFreq]
        ]
      centerR2Freq = div numR2Freqs 2
      radialConstant = 2 * pi / (log periodEnv)
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             VG.convert . toUnboxed . computeS $
             fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. x :. y) ->
               let xFreq = fromIntegral $ x - centerR2Freq
                   yFreq = fromIntegral $ y - centerR2Freq
                   rho = sqrt $ xFreq ^ 2 + yFreq ^ 2
                   phi = atan2 yFreq xFreq
               in if rho == 0
                    then 0
                    else ((0 :+ (-1)) ^ 1) *
                         ((periodR2 / (pi * rho) :+ 0) **
                          (0 :+ fromIntegral radialFreq * radialConstant)) *
                         (cis $ phi * fromIntegral angularFreq))
          idxs
  in pinwheels

  
{-# INLINE analyticalFourierSeries2' #-}
analyticalFourierSeries2' ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierSeries2' numPoints _ angularFreq radialFreq sigma periodR2 periodEnv =
  let center = div numPoints 2
      arr =
        fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
          let x' = i - center
              y' = j - center
          in if x' == 0 && y' == 0
               then 0
               else let x =
                          if x' == 0
                            then 0
                            else if x' > 0
                                   then 1 / (fromIntegral $ center + 1 - x')
                                   else (-1) / (fromIntegral $ center + 1 + x')
                        y =
                          if y' == 0
                            then 0
                            else if y' > 0
                                   then 1 / (fromIntegral $ center + 1 - y')
                                   else (-1) / (fromIntegral $ center + 1 + y')
                        rho = sqrt $ x ^ 2 + y ^ 2
                        phi = atan2 y x
                    in analyticalFourierSeriesFunc2
                         angularFreq
                         radialFreq
                         sigma
                         periodR2
                         periodEnv
                         phi
                         rho
  in arr
  
{-# INLINE analyticalFourierCoefficients2' #-}
analyticalFourierCoefficients2' ::
     (Eq e, Fractional e, RealFloat e, Gamma (Complex e), Unbox e)
  => Int
  -> e
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> R.Array D DIM2 (Complex e)
analyticalFourierCoefficients2' numFreqs delta angularFreq radialFreq sigma periodR2 periodEnv =
  let c = (-1) ^ (abs angularFreq) :+ 0
  in R.map (* c) $
     analyticalFourierSeries2'
       numFreqs
       delta
       angularFreq
       radialFreq
       sigma
       periodR2
       periodEnv


pinwheelFourierCoefficientsAnatical' ::
     ( Unbox e
     , RealFloat e
     , Gamma (Complex e)
     , VG.Vector vector (Complex e)
     , NFData (vector (Complex e))
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> e
  -> e
  -> IA.Array (Int, Int) (vector (Complex e))
pinwheelFourierCoefficientsAnatical' numR2Freqs phiFreq rhoFreq thetaFreq rFreq sigma periodR2 periodEnv =
  let idxs =
        [ (radialFreq, angularFreq)
        | radialFreq <- pinwheelFreqs rhoFreq rFreq
        , angularFreq <- pinwheelFreqs phiFreq thetaFreq
        ]
      centerR2Freq = div numR2Freqs 2
      (radialLB, radialUB) = pinwheelFreqsBound rhoFreq rFreq
      (angularLB, angularUB) = pinwheelFreqsBound phiFreq thetaFreq
      pinwheels =
        parMap
          rdeepseq
          (\(radialFreq, angularFreq) ->
             VG.convert .
             toUnboxed .
             computeS -- .
             -- R.zipWith
             --   (*)
             --   (fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
             --      gaussian2D
             --        (fromIntegral $ i - centerR2Freq)
             --        (fromIntegral $ j - centerR2Freq)
             --        (fromIntegral $ div centerR2Freq 2))
            $
             analyticalFourierCoefficients2'
               numR2Freqs
               1
               (angularFreq)
               (radialFreq)
               sigma
               periodR2
               periodEnv)
          idxs
  in listArray ((radialLB, angularLB), (radialUB, angularUB)) pinwheels
  
{-# INLINE idealHighPassFilter #-}
idealHighPassFilter ::
     (VG.Vector vector (Complex Double))
  => Double
  -> Int
  -> vector (Complex Double)
idealHighPassFilter radius numR2Freq =
  let a = 1 / (radius * 2)
      r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      sinc x =
        if x == 0
          then 1
          else (sin (pi * x)) / (pi * x)
  in VG.fromList
       [ if xFreq == 0 && yFreq == 0
         then (1 - 1 / a ^ 2) :+ 0
         else ((-1) * (sinc (xFreq / a)) * (sinc (yFreq / a)) / (a ^ 2)) :+ 0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]

{-# INLINE idealHighPassFilter1 #-}
idealHighPassFilter1 :: Int -> VS.Vector (Complex Double)
idealHighPassFilter1 numR2Freq =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
  in VS.fromList
       [ if xFreq == 0 && yFreq == 0
         then 0
         else (-1) :+ 0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]
  
{-# INLINE gaussianHighPassFilter #-}
gaussianHighPassFilter :: Double -> Int -> VS.Vector (Complex Double)
gaussianHighPassFilter alpha numR2Freq =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
  in VS.fromList
       [ if xFreq == 0 && yFreq == 0
         then (1 - (pi / alpha)) :+ 0
         else ((-1) * (pi / alpha) *
               exp ((pi ^ 2 * (xFreq ^ 2 + yFreq ^ 2)) / (-alpha))) :+
              0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]
  
{-# INLINE gaussianHighPassFilter1 #-}
gaussianHighPassFilter1 :: Double -> Double -> Int -> VS.Vector (Complex Double)
gaussianHighPassFilter1 radius alpha numR2Freq =
  let a = 1 / (radius * 2)
      r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      sinc x =
        if x == 0
          then 1
          else (sin (pi * x)) / (pi * x)
  in VS.fromList
       [ (sinc (xFreq / a) * sinc (yFreq / a) / (a ^ 2) -
          (pi / alpha) * exp ((pi ^ 2 * (xFreq ^ 2 + yFreq ^ 2)) / (-alpha))) :+
       0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]
  
{-# INLINE idealLowPassFilter #-}
idealLowPassFilter ::
     (VG.Vector vector (Complex Double))
  => Double
  -> Double
  -> Int
  -> vector (Complex Double)
idealLowPassFilter radius periodR2 numR2Freq =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
   in VG.fromList
        [ let rho = 2 * pi * sqrt (xFreq ^ 2 + yFreq ^ 2) / periodR2
           in if rho == 0
                then 0
                else radius / rho  * bessel_J1 (radius * rho) / periodR2 / (2 * pi)^2 :+ 0
        | xFreq <- r2Freqs
        , yFreq <- r2Freqs
        ]

{-# INLINE gaussianLowPassFilter #-}
gaussianLowPassFilter :: Double -> Int -> VS.Vector (Complex Double)
gaussianLowPassFilter std numR2Freq =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
  in VS.fromList
       [ ((std ^ 2) * exp ((-pi) * (xFreq ^ 2 + yFreq ^ 2) * (std ^ 2))) :+
       0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]

{-# INLINE laplacianLowPassFilter #-}
laplacianLowPassFilter :: Double -> Int -> VS.Vector (Complex Double)
laplacianLowPassFilter a numR2Freq =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
  in VS.fromList
       [ let s2 = xFreq ^ 2 + yFreq ^ 2
         in if xFreq == 0 && yFreq == 0
              then (1 - a / ((4 * pi * pi * s2 + a ^ 2) ** 1.5)) :+ 0
              else (-a) / ((4 * pi * pi * s2 + a ^ 2) ** 1.5) :+ 0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]
       
{-# INLINE laplacianHighPassFilter #-}
laplacianHighPassFilter :: Int -> VS.Vector (Complex Double)
laplacianHighPassFilter numR2Freq =
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
  in VS.fromList
       [ if xFreq == 0 && yFreq == 0
         then 8 -- -3.33
         else if xFreq == 0 && abs yFreq == 1
                then -1 -- 0.67
                else if abs xFreq == 1 && yFreq == 0
                       then -1 -- 0.67
                       else if abs xFreq == 1 && abs yFreq == 1
                              then -1 -- 0.17
                              else 0
       | xFreq <- r2Freqs
       , yFreq <- r2Freqs
       ]

{-# INLINE convolveFrequency #-}
convolveFrequency ::
     (VG.Vector vector (Complex Double))
  => DFTPlan
  -> Int
  -> vector (Complex Double)
  -> IA.Array (Int, Int) (vector (Complex Double))
  -> IO (IA.Array (Int, Int) (vector (Complex Double)))
convolveFrequency plan numR2Freq filter' arr = do
  let dftPlanID = DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]
      idftPlanID = DFTPlanID IDFT1DG [numR2Freq, numR2Freq] [0, 1]
      vecs = L.map VG.convert . IA.elems $ arr
      filter =
        VU.convert .
        toUnboxed .
        computeUnboxedS .
        makeFilter2D . fromUnboxed (Z :. numR2Freq :. numR2Freq) . VG.convert $
        filter'
  dftF <- dftExecute plan dftPlanID filter
  vecsF <- mapConcurrently (dftExecute plan dftPlanID) vecs
  outputs <-
    mapConcurrently (dftExecute plan idftPlanID . VS.zipWith (*) dftF) vecsF
  return . listArray (bounds arr) . L.map VG.convert $ outputs
  
{-# INLINE convolveFrequency1 #-}
convolveFrequency1 ::
     (VG.Vector vector (Complex Double))
  => DFTPlan
  -> Int
  -> vector (Complex Double)
  -> IA.Array (Int, Int) (vector (Complex Double))
  -> IO (IA.Array (Int, Int) (vector (Complex Double)))
convolveFrequency1 plan numR2Freq filter' arr = do
  let dftPlanID = DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]
      idftPlanID = DFTPlanID IDFT1DG [numR2Freq, numR2Freq] [0, 1]
      vecs = L.map VG.convert . IA.elems $ arr
      filter =
        VU.convert .
        toUnboxed .
        computeUnboxedS .
        makeFilter2D . fromUnboxed (Z :. numR2Freq :. numR2Freq) . VG.convert $
        filter'
  dftF <- dftExecute plan dftPlanID filter
  vecsF <- mapConcurrently (dftExecute plan dftPlanID) vecs
  outputs <-
    mapConcurrently (dftExecute plan idftPlanID . VS.zipWith (*) dftF) vecsF
  return .
    listArray (bounds arr) . L.map VG.convert . L.zipWith (VG.zipWith (-)) vecs $
    outputs
    
{-# INLINE centerHollow #-}
centerHollow ::
     (VG.Vector vector (Complex Double))
  => Int
  -> IA.Array (Int, Int) (vector (Complex Double))
  -> IA.Array (Int, Int) (vector (Complex Double))
centerHollow numR2Freq arr = centerHollowVector numR2Freq <$> arr

{-# INLINE centerHollowVector #-}
centerHollowVector ::
     (VG.Vector vector (Complex Double))
  => Int
  -> vector (Complex Double)
  -> vector (Complex Double)
centerHollowVector numR2Freq vec =
  let s = VG.sum vec / (fromIntegral (numR2Freq ^ 2) :+ 0)
   in VG.map (\x -> x - s) vec   

{-# INLINE centerHollowArray #-}
centerHollowArray ::
     (R.Source s (Complex e), Unbox e, RealFloat e)
  => Int
  -> R.Array s DIM2 (Complex e)
  -> R.Array D DIM2 (Complex e)
centerHollowArray numR2Freq arr =
  let s = sumAllS arr / (fromIntegral (numR2Freq ^ 2) :+ 0)
   in R.map (\x -> x - s) $ arr
   
{-# INLINE centerHollowArray' #-}
centerHollowArray' ::
     (R.Source s (Complex e), Unbox e, RealFloat e)
  => Int
  -> R.Array s DIM2 (Complex e)
  -> R.Array D DIM2 (Complex e)
centerHollowArray' numR2Freq arr =
  let s = sumAllS arr
      c = div numR2Freq 2
   in R.traverse arr id $ \f idx@(Z :. i :. j) ->
        if i == c && j == c
          then (-s)
          else f idx

envelopIntegral ::
     Double -> Double -> Double -> Double -> Double -> Int -> Complex Double
envelopIntegral a b delta s period radialFreq =
  let m = round $ (b - a) / delta
      n =
        if odd m
          then m
          else m - 1
      -- weights = VU.fromList $ weightsSimpsonRule n
      vec =
        VU.generate n $ \i ->
          let x = (a + fromIntegral i * delta)
          -- in   ((x :+ 0) ** (s :+ freq)) * ((pi :+ 0) ** (0 :+ freq)) *
          --      (gamma ((0.5 :+ 0) * (1 :+ (-freq)))) /
          --      (gamma ((0.5 :+ 0) * (1 :+ freq)))
          in analyticalFourierSeriesFunc2 0 radialFreq s period period 0 x
  -- in (delta / 3 :+ 0) * (VU.sum . VU.zipWith (*) vec $ weights)
  in VU.sum vec

envelopIntegral2D ::
     Int -> Double -> Double -> Double -> Int -> Complex Double
envelopIntegral2D numPoints delta s period radialFreq =
  let center = div numPoints 2
      -- weights = VU.fromList $ weightsSimpsonRule n
      arr = analyticalFourierSeries2 numPoints delta 0 radialFreq s period period
  in sumAllS arr

printEnvelopIntegral ::
     Double -> Double -> Double -> Double -> Double -> Int -> IO ()
printEnvelopIntegral a b delta s period radialFreq =
  let m = round $ (b - a) / delta
      n =
        if odd m
          then m
          else m - 1
      weights = VU.fromList $ weightsSimpsonRule n
      vec =
        VU.generate n $ \i ->
          let x = a + fromIntegral i * delta  
          in analyticalFourierSeriesFunc2 0 radialFreq s period period 0 x
  in print $ VU.zipWith (*) weights vec
