{-# LANGUAGE FlexibleContexts #-}
module STC.Plan where

import           Control.Monad        as M
import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           DFT.Plan
import           System.Random
import           Types

makeR2Z1T0Plan ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> R.Array r DIM4 (Complex Double)
  -> IO DFTPlan
makeR2Z1T0Plan oldPlan arr = do
  let (Z :. numThetaFreqs :. numTheta0Freqs :. xLen :. yLen) = extent arr
  lock <- getFFTWLock
  vecTemp1 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp2 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (numThetaFreqs * xLen * yLen) randomIO :: IO (VS.Vector (Complex Double))
  vecTemp3 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (xLen * yLen * numTheta0Freqs) randomIO :: IO (VS.Vector (Complex Double))
  fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [numThetaFreqs, numTheta0Freqs, xLen, yLen]
       [2, 3]
       vecTemp1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [numThetaFreqs, numTheta0Freqs, xLen, yLen]
         [2, 3]
         vec >>= \(plan, _) ->
         dft1dGPlan lock plan [numThetaFreqs, xLen, yLen] [0] vecTemp2 >>= \(plan, vec) ->
           idft1dGPlan lock plan [numThetaFreqs, xLen, yLen] [0] vec >>= \(plan, _) ->
             dft1dGPlan lock plan [numTheta0Freqs, xLen, yLen] [1, 2] vecTemp3 >>= \(plan, vec) ->
               idft1dGPlan lock plan [numTheta0Freqs, xLen, yLen] [1, 2] vec)


makeR2Z2T0S0Plan ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> Bool
  -> FilePath
  -> R.Array r DIM6 (Complex Double)
  -> IO DFTPlan
makeR2Z2T0S0Plan oldPlan wisdomFlag wisdomFilePath arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. numTheta0Freqs :. numScale0Freqs :. xLen :. yLen) =
        extent arr
  when wisdomFlag (importFFTWWisdom wisdomFilePath)
  lock <- getFFTWLock
  let vecTemp1 = VS.convert . toUnboxed . computeS . delay $ arr
  plan <-
    fst <$>
    (-- dft1dGPlan
     --   lock
     --   oldPlan
     --   [ numThetaFreqs
     --   , numScaleFreqs
     --   , numTheta0Freqs
     --   , numScale0Freqs
     --   , xLen
     --   , yLen
     --   ]
     --   [4, 5]
     --   vecTemp1 >>= \(plan, vec) ->
     --   idft1dGPlan
     --     lock
     --     plan
     --     [ numThetaFreqs
     --     , numScaleFreqs
     --     , numTheta0Freqs
     --     , numScale0Freqs
     --     , xLen
     --     , yLen
     --     ]
     --     [4, 5]
     --     vec >>= \(plan, vec) ->
         dft1dGPlan
           lock
           oldPlan
           [numThetaFreqs, numScaleFreqs, xLen, yLen]
           [0, 1]
           vecTemp1 >>= \(plan, vec) ->
           idft1dGPlan
             lock
             plan
             [numThetaFreqs, numScaleFreqs, xLen, yLen]
             [0, 1]
             vec >>= \(plan, vec) ->
             dft1dGPlan
               lock
               plan
               [numTheta0Freqs, numScale0Freqs, xLen, yLen]
               [2, 3]
               vec >>= \(plan, vec) ->
               idft1dGPlan
                 lock
                 plan
                 [numTheta0Freqs, numScale0Freqs, xLen, yLen]
                 [2, 3]
                 vec >>= \(plan, vec) ->
                 dft1dGPlan lock plan [numScale0Freqs, xLen, yLen] [1, 2] vec)
  exportFFTWWisdom wisdomFilePath
  return plan


makeImagePlan ::
     DFTPlan
  -> R.Array U DIM3 Double
  -> IO (DFTPlan, R.Array U DIM3 (Complex Double))
makeImagePlan plan arr = do
  let (Z :. channels :. cols :. rows) = extent arr
  lock <- getFFTWLock
  (plan1, imgF) <-
    dft1dGPlan lock plan [channels, cols, rows] [1, 2] .
    VS.map (:+ 0) . VS.convert . toUnboxed $
    arr
  (newPlan, _) <- idft1dGPlan lock plan1 [channels, cols, rows] [1, 2] imgF
  return (newPlan, fromUnboxed (extent arr) . VS.convert $ imgF)


-- For local eigenvetor

make4DPlan ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> Bool
  -> FilePath
  -> R.Array r DIM4 (Complex Double)
  -> IO DFTPlan
make4DPlan oldPlan wisdomFlag wisdomFilePath arr = do
  let (Z :. numThetaFreqs :. numScaleFreqs :. xLen :. yLen) = extent arr
  when wisdomFlag (importFFTWWisdom wisdomFilePath)
  lock <- getFFTWLock
  let vecTemp1 = VS.convert . toUnboxed . computeS . delay $ arr 
  plan <-
    fst <$>
    (dft1dGPlan
       lock
       oldPlan
       [numThetaFreqs, numScaleFreqs, xLen, yLen]
       [2,3]
       vecTemp1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [numThetaFreqs, numScaleFreqs, xLen, yLen]
         [2,3]
         vec >>= \(plan, vec) ->
         dft1dGPlan lock plan [xLen, yLen] [0, 1] vec >>= \(plan, vec) ->
           dft1dGPlan
             lock
             plan
             [numThetaFreqs, numScaleFreqs, xLen, yLen]
             [0, 1]
             vecTemp1 >>= \(plan, vec) ->
             idft1dGPlan
               lock
               plan
               [numThetaFreqs, numScaleFreqs, xLen, yLen]
               [0, 1]
               vec)
  exportFFTWWisdom wisdomFilePath
  return plan
