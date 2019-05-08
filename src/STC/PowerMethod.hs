{-# LANGUAGE FlexibleContexts #-}
module STC.PowerMethod
  ( module STC.PowerMethod
  , module STC.InitialDistribution
  , module STC.Plan
  , module STC.Bias
  , module STC.Utils
  ) where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           FokkerPlanck.DomainChange (r2s1tor2z1, r2z1Tor2s1,
                                            r2z2Tor2s1rp, r2z2Tor2s1rpP)
import           STC.Convolution
import           Types
import           Utils.Array

import           Filter.Utils
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           STC.CompletionField
import           STC.InitialDistribution
import STC.Plan
import           System.FilePath           ((</>))
import           Text.Printf
import           Utils.Time
import STC.Bias
import STC.Utils
import STC.Reversal

{-# INLINE eigenVectorR2Z1 #-}
eigenVectorR2Z1 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> R2Z1T0Array
  -> Int
  -> Bool
  -> String
  -> R.Array s1 DIM3 (Complex Double)
  -> R.Array s2 DIM4 (Complex Double)
  -> IO R2Z1T0Array
eigenVectorR2Z1 plan folderPath numOrientation thetaFreqs filterF n writeFlag name bias inputR2Z1T0 = do
  let (Z :. numThetaFreq :. numTheta0Freq :. cols :. rows) = extent inputR2Z1T0
  inputR2Z1 <- R.sumP . rotateR2Z1T0Array $ inputR2Z1T0
  let sourceDist = R.zipWith (*) bias inputR2Z1
      s = VU.maximum . VU.map magnitude . toUnboxed . computeS $ sourceDist
  let initialDist = computeS $ R.map (\x -> x / (s :+ 0)) sourceDist
  sourceArr <- convolveR2T0 plan filterF initialDist
  when
    writeFlag
    (do sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceArr
        sourceField <-
          computeP .
          R.extend (Z :. (1 :: Int) :. All :. All) .
          R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
          sourceR2Z1
        plotImageRepaComplex
          (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
          ImageRepa 8 $
          sourceField)
  return sourceArr

-- eigensink is computed from eigensource
powerMethod1 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> R2Z1T0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array s1 DIM3 (Complex Double)
  -> R.Array s2 DIM4 (Complex Double)
  -> IO (R.Array D DIM3 (Complex Double))
powerMethod1 plan folderPath cols rows numOrientation thetaFreqs theta0Freqs filter numIteration writeFlag idStr threshold bias eigenVecSource = do
  filterF <- dftR2Z1T0 plan . computeS . makeFilter2D $ filter
  sourceR2Z1T0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z1
           plan
           folderPath
           numOrientation
           thetaFreqs
           filterF
           n
           writeFlag
           "Source"
           bias
           input)
      (computeS . delay $ eigenVecSource)
      [1 .. numIteration]
  let sinkR2Z1T0 =
        computeSinkFromSourceR2Z1T0 thetaFreqs theta0Freqs sourceR2Z1T0
  sourceR2Z1 <- R.sumP . rotateR2Z1T0Array $ sourceR2Z1T0
  sinkR2Z1 <- R.sumP . rotateR2Z1T0Array $ sinkR2Z1T0
  sinkField <-
    computeP .
    R.extend (Z :. (1 :: Int) :. All :. All) .
    R.sumS . rotate3D . r2z1Tor2s1 numOrientation thetaFreqs $
    sinkR2Z1
  plotImageRepaComplex (folderPath </> "Sink.png") . ImageRepa 8 $ sinkField
  completionFieldR2 <-
    completionFieldR2Z1
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      sourceR2Z1
      sinkR2Z1
  let completionFieldR2' =
        R.zipWith (\x y -> x * magnitude y) completionFieldR2 . R.slice bias $
        (Z :. (0 :: Int) :. All :. All)
      m = threshold * (VU.maximum . toUnboxed . computeS $ completionFieldR2)
      newBias =
        R.traverse2 bias completionFieldR2 const $ \fb fc idx@(Z :. _ :. i :. j) ->
          if fc (Z :. i :. j) >= m
            then 0
            else fb idx
  -- plotImageRepa
  --   (folderPath </> printf "bias%s.png" idStr)
  --   (ImageRepa 8 .
  --    computeS .
  --    R.map magnitude .
  --    R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
  --    newBias)
  return newBias

{-# INLINE eigenVectorR2Z2 #-}
eigenVectorR2Z2 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2 plan folderPath numOrientation thetaFreqs numScale scaleFreqs filterF n writeFlag name bias inputR2Z2T0S0 = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
        extent inputR2Z2T0S0
  inputR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ inputR2Z2T0S0) >>= R.sumP
  let sourceDist = R.zipWith (*) bias inputR2Z2
  s <- R.foldAllP max 0 . R.map magnitude $ sourceDist
  let initialDist = R.map (/ (s :+ 0)) sourceDist
  sourceArr <- convolveR2T0S0P plan filterF initialDist
  printCurrentTime $ printf "iteration %d" n
  when
    (n == 1 || (writeFlag && odd n))
    (let
      in do sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceArr) >>= R.sumP
            sourceR2S1RP <-
              r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
              sourceR2Z2
            sourceField <-
              (R.sumP . rotate4D . rotate4D $ sourceR2S1RP) >>=
              fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
              R.sumP
            plotImageRepaComplex
              (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
              ImageRepa 8 $
              sourceField)
  return sourceArr

powerMethodR2Z2T0S0 ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> Int
  -> [Double]
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0 plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs filter numIteration writeFlag idStr threshold bias eigenVecSource = do
  filterF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filter
  sourceR2Z2T0S0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           filterF
           n
           writeFlag
           ("Source" L.++ idStr)
           bias
           input)
      (computeS . delay $ eigenVecSource)
      [1 .. numIteration]
  let sinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceR2Z2T0S0
  sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceR2Z2T0S0) >>= R.sumP
  sinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sinkR2Z2T0S0) >>= R.sumP
  sinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     sinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s.png" idStr) . ImageRepa 8 $
    sinkField
  completionFieldR2 <-
    completionFieldR2Z2
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      sourceR2Z2
      sinkR2Z2
  let m = threshold * (VU.maximum . toUnboxed . computeS $ completionFieldR2)
      newBias =
        R.traverse2 bias completionFieldR2 const $ \fb fc idx@(Z :. _ :. _ :. i :. j) ->
          if fc (Z :. i :. j) >= m
            then 0
            else fb idx
  -- plotImageRepa
  --   (folderPath </> printf "bias%s.png" idStr)
  --   (ImageRepa 8 .
  --    computeS .
  --    R.map magnitude .
  --    R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
  --    newBias)
  return newBias    



{-# INLINE eigenVectorR2Z2Bias #-}
eigenVectorR2Z2Bias ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2Bias plan folderPath numOrientation thetaFreqs numScale scaleFreqs filterF n writeFlag name bias inputR2Z2T0S0 = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
        extent inputR2Z2T0S0
  inputR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ inputR2Z2T0S0) >>= R.sumP
  sourceDist <- convolveR2Z2 plan bias inputR2Z2
  s <- R.foldAllP max 0 . R.map magnitude $ sourceDist
  let initialDist = R.map (/ (s :+ 0)) sourceDist
  sourceArr <- convolveR2T0S0P plan filterF initialDist
  printCurrentTime $ printf "iteration %d" n
  when
    (n == 1 || (writeFlag && odd n))
    (let
      in do sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceArr) >>= R.sumP
            sourceR2S1RP <-
              r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
              sourceR2Z2
            sourceField <-
              (R.sumP . rotate4D . rotate4D $ sourceR2S1RP) >>=
              fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
              R.sumP
            plotImageRepaComplex
              (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
              ImageRepa 8 $
              sourceField)
  return sourceArr


powerMethodR2Z2T0S0Bias ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> Int
  -> [Double]
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0Bias plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs filter numIteration writeFlag idStr threshold bias eigenVecSource = do
  filterF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filter
  sourceR2Z2T0S0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2Bias
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           filterF
           n
           writeFlag
           ("Source" L.++ idStr)
           bias
           input)
      (computeS $ delay eigenVecSource)
      [1 .. numIteration]
  let sinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceR2Z2T0S0
  sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceR2Z2T0S0) >>= R.sumP
  sinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sinkR2Z2T0S0) >>= R.sumP
  sinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     sinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s.png" idStr) . ImageRepa 8 $
    sinkField
  completionFieldR2 <-
    completionFieldR2Z2
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      sourceR2Z2
      sinkR2Z2
  let m = threshold * (VU.maximum . toUnboxed . computeS $ completionFieldR2)
      newBias =
        R.traverse2 bias completionFieldR2 const $ \fb fc idx@(Z :. _ :. _ :. i :. j) ->
          if fc (Z :. i :. j) >= m
            then 0
            else fb idx
  -- plotImageRepa
  --   (folderPath </> printf "bias%s.png" idStr)
  --   (ImageRepa 8 .
  --    computeS .
  --    R.map magnitude .
  --    R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
  --    newBias)
  return newBias


{-# INLINE eigenVectorR2Z2Reversal #-}
eigenVectorR2Z2Reversal ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2Reversal plan folderPath numOrientation thetaFreqs numScale scaleFreqs filterF n writeFlag name reversalFactor bias inputR2Z2T0S0' = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
        extent inputR2Z2T0S0'
      reversal = computeReversalR2Z2T0S0 thetaFreqs inputR2Z2T0S0'
      inputR2Z2T0S0 =
        R.zipWith
          (\x y -> x + (reversalFactor :+ 0) * y)
          inputR2Z2T0S0'
          reversal
  inputR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ inputR2Z2T0S0) >>= R.sumP
  let sourceDist = R.zipWith (*) bias inputR2Z2
  s <- R.foldAllP max 0 . R.map magnitude $ sourceDist
  let initialDist = R.map (/ (s :+ 0)) sourceDist
  sourceArr <- convolveR2T0S0P plan filterF initialDist
  printCurrentTime $ printf "iteration %d" n
  when
    (n == 1 || (writeFlag && odd n))
    (let
      in do sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceArr) >>= R.sumP
            sourceR2S1RP <-
              r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
              sourceR2Z2
            sourceField <-
              (R.sumP . rotate4D . rotate4D $ sourceR2S1RP) >>=
              fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
              R.sumP
            plotImageRepaComplex
              (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
              ImageRepa 8 $
              sourceField)
  return sourceArr

powerMethodR2Z2T0S0Reversal ::
     (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> [Double]
  -> [Double]
  -> Int
  -> [Double]
  -> [Double]
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> Double
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0Reversal plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs filter numIteration writeFlag idStr threshold reversalFactor bias eigenVecSource = do
  filterF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filter
  sourceR2Z2T0S0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2Reversal
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           filterF
           n
           writeFlag
           ("Source" L.++ idStr)
           reversalFactor
           bias
           input)
      (computeS . delay $ eigenVecSource)
      [1 .. numIteration]
  let sinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceR2Z2T0S0
  sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceR2Z2T0S0) >>= R.sumP
  sinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sinkR2Z2T0S0) >>= R.sumP
  sinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     sinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s.png" idStr) . ImageRepa 8 $
    sinkField
  completionFieldR2 <-
    completionFieldR2Z2
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      sourceR2Z2
      sinkR2Z2
  let m = threshold * (VU.maximum . toUnboxed . computeS $ completionFieldR2)
      newBias =
        R.traverse2 bias completionFieldR2 const $ \fb fc idx@(Z :. _ :. _ :. i :. j) ->
          if fc (Z :. i :. j) >= m
            then 0
            else fb idx
  -- plotImageRepa
  --   (folderPath </> printf "bias%s.png" idStr)
  --   (ImageRepa 8 .
  --    computeS .
  --    R.map magnitude .
  --    R.extend (Z :. (1 :: Int) :. All :. All) . R.sumS . rotate3D $
  --    newBias)
  return newBias
