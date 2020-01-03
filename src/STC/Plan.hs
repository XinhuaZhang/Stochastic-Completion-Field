{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
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
import           System.Random

{-# INLINE makePlan #-}
makePlan :: DFTPlan -> Int -> Int -> Int -> Int -> IO DFTPlan
makePlan !initPlan !nx !ny !numThetaFreq !numRFreq = do
  initVec1 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$> M.replicateM (nx * ny) randomIO
  initVec2 <-
    (VS.fromList . L.map (\x -> x :+ 0)) <$>
    M.replicateM (nx * ny * numThetaFreq * numRFreq) randomIO
  lock <- getFFTWLock
  fst <$>
    (dft1dGPlan lock initPlan [nx, ny] [0, 1] initVec1 >>= \(plan, vec) ->
       idft1dGPlan lock plan [nx, ny] [0, 1] vec >>= \(plan, vec) ->
         dft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] initVec2 >>= \(plan, vec) ->
           idft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec)

{-# INLINE makePlanFromArray #-}
makePlanFromArray ::
     (R.Source r (Complex Double))
  => DFTPlan
  -> R.Array r DIM4 (Complex Double)
  -> IO DFTPlan
makePlanFromArray !initPlan arr = do
  let (Z :. (!numRFreq) :. (!numThetaFreq) :. (!nx) :. (!ny)) = extent arr
      !initVec = VS.fromList . R.toList $ arr
  lock <- getFFTWLock
  fst <$>
    (dft1dGPlan lock initPlan [nx, ny] [0, 1] initVec >>= \(plan, vec) ->
       idft1dGPlan lock plan [nx, ny] [0, 1] vec >>= \(plan, vec) ->
         dft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec >>= \(plan, vec) ->
           idft1dGPlan lock plan [numRFreq, numThetaFreq, nx, ny] [0, 1] vec)
