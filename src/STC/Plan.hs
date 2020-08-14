{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict #-}
module STC.Plan
  ( module DFT.Plan
  , module STC.Plan
  ) where

import           Control.Monad        as M
import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           DFT.Plan
import           System.FilePath
import           System.Random

{-# INLINE makePlan #-}
makePlan :: FilePath -> DFTPlan -> Int -> Int -> Int -> Int -> Int -> Int -> Int -> IO DFTPlan
makePlan folderPath initPlan nx ny numR2Freq numThetaFreq numRFreq numPhiFreq numRhoFreq = do
  importFFTWWisdom (folderPath </> "fftwwisdom.dat")
  initVec1 <-
    VS.fromList . L.map (:+ 0) <$>
    M.replicateM (numRFreq * numThetaFreq * numR2Freq ^ 2) randomIO
  initVec2 <-
    VS.fromList . L.map (:+ 0) <$>
    M.replicateM (nx * ny * numThetaFreq * numRFreq) randomIO
  initVec3 <-
    VS.fromList . L.map (:+ 0) <$> M.replicateM (numR2Freq ^ 2) randomIO
  initVec4 <-
    VS.fromList . L.map (:+ 0) <$>
    M.replicateM (numThetaFreq * numRFreq * numPhiFreq * numRhoFreq) randomIO
  initVec5 <-
    VS.fromList . L.map (:+ 0) <$> M.replicateM (numRhoFreq * numRFreq) randomIO
  lock <- getFFTWLock
  plan <-
    fst <$>
    (dft1dGPlan
       lock
       initPlan
       [numRFreq, numThetaFreq, numR2Freq, numR2Freq]
       [0, 1, 2, 3]
       initVec1 >>= \(plan, vec) ->
       idft1dGPlan
         lock
         plan
         [numRFreq, numThetaFreq, numR2Freq, numR2Freq]
         [0, 1, 2, 3]
         vec >>= \(plan, _) ->
         dft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] initVec2 >>= \(plan, vec) ->
           idft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec >>= \(plan, _) ->
             dft1dGPlan lock plan [numR2Freq, numR2Freq] [0, 1] initVec3 >>= \(plan, vec) ->
               idft1dGPlan lock plan [numR2Freq, numR2Freq] [0, 1] vec >>= \(plan, _) ->
                 dft1dGPlan
                   lock
                   plan
                   [numRFreq, numRhoFreq, numThetaFreq, numPhiFreq]
                   [0, 1]
                   initVec4 >>= \(plan, vec) ->
                   idft1dGPlan
                     lock
                     plan
                     [numRFreq, numRhoFreq, numThetaFreq, numPhiFreq]
                     [0, 1]
                     vec >>= \(plan, _) ->
                     dft1dGPlan lock plan [numRFreq, numRhoFreq] [0, 1] initVec5 >>= \(plan, vec) ->
                       idft1dGPlan lock plan [numRFreq, numRhoFreq] [0, 1] vec)
  exportFFTWWisdom (folderPath </> "fftwwisdom.dat")
  return plan

{-# INLINE makePlanFromArray #-}
makePlanFromArray ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> R.Array r DIM4 (Complex Double)
  -> IO DFTPlan
makePlanFromArray initPlan arr = do
  let (Z :. (numRFreq) :. (numThetaFreq) :. (nx) :. (ny)) = extent arr
      initVec = VS.fromList . R.toList $ arr
  lock <- getFFTWLock
  fst <$>
    (dft1dGPlan lock initPlan [nx, ny] [0, 1] initVec >>= \(plan, vec) ->
       idft1dGPlan lock plan [nx, ny] [0, 1] vec >>= \(plan, vec) ->
         dft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec >>= \(plan, vec) ->
           idft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec)


-- {-# INLINE makePlan #-}
-- makePlan :: FilePath -> DFTPlan -> Int -> Int -> Int -> Int -> Int -> IO DFTPlan
-- makePlan folderPath initPlan nx ny numR2Freq numThetaFreq numRFreq = do
--   importFFTWWisdom (folderPath </> "fftwwisdom.dat")
--   initVec1 <-
--     VS.fromList . L.map (:+ 0) <$> M.replicateM (numR2Freq ^ 2) randomIO
--   initVec2 <-
--     VS.fromList . L.map (:+ 0) <$>
--     M.replicateM (nx * ny * numThetaFreq * numRFreq) randomIO
--   initVec3 <- VS.fromList . L.map (:+ 0) <$> M.replicateM (nx * ny) randomIO
--   lock <- getFFTWLock
--   plan <-
--     fst <$>
--     (dft1dGPlan lock initPlan [numR2Freq, numR2Freq] [0, 1] initVec1 >>= \(plan, vec) ->
--        idft1dGPlan lock plan [numR2Freq, numR2Freq] [0, 1] vec >>= \(plan, _) ->
--          dft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] initVec2 >>= \(plan, vec) ->
--            idft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec >>= \(plan, _) ->
--              dft1dGPlan lock plan [nx, ny] [0, 1] initVec3 >>= \(plan, vec) ->
--                idft1dGPlan lock plan [nx, ny] [0, 1] vec)
--   exportFFTWWisdom (folderPath </> "fftwwisdom.dat")
--   return plan
