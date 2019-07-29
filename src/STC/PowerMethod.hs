{-# LANGUAGE FlexibleContexts #-}
module STC.PowerMethod where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Complex
import           Data.List                 as L
import           Data.Vector.Storable      as VS
import           Data.Vector.Unboxed       as VU
import           DFT.Plan
import           FokkerPlanck.DomainChange (r2s1tor2z1, r2z1Tor2s1,
                                            r2z2Tor2s1rp, r2z2Tor2s1rpP)
import           STC.Convolution
import           Types
import           Utils.Array

import           Data.Ix
import           Filter.Utils
import           Image.IO                  (ImageRepa (..), plotImageRepa,
                                            plotImageRepaComplex)
import           Image.Transform
import           STC.Bias
import           STC.CompletionField
import           STC.InitialDistribution
import           STC.Plan
import           STC.PowerMethodNormalization
import           STC.Reversal
import           STC.Utils
import           System.FilePath           ((</>))
import           Text.Printf
import           Utils.Time

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
  -> Double
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2 plan folderPath numOrientation thetaFreqs numScale scaleFreqs maxScale filterF n writeFlag name bias inputR2Z2T0S0 = do
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
              r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs maxScale $
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
  -> Double
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array s1 DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
powerMethodR2Z2T0S0 plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs maxScale filter numIteration writeFlag idStr threshold bias eigenVecSource = do
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
           maxScale
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

{-# INLINE eigenVectorR2Z2Bias #-}
eigenVectorR2Z2Bias ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> Double
  -> R2Z2T0S0Array
  -> PowerMethodNormalizationOption
  -> Int
  -> Bool
  -> String
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2Bias plan folderPath numOrientation thetaFreqs numScale scaleFreqs maxScale filterF normMethod n writeFlag name bias inputR2Z2T0S0 = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
        extent inputR2Z2T0S0
  inputR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ inputR2Z2T0S0) >>= R.sumP
  biasedInputR2Z2 <- convolveR2Z2 plan bias inputR2Z2
  normalizedBiasedInputR2Z2 <-
    powerMethodNormalization normMethod biasedInputR2Z2
  sourceArr <- convolveR2T0S0P plan filterF normalizedBiasedInputR2Z2
  printCurrentTime $ printf "iteration %d" (n + 1)
  when
    (n == 0 || writeFlag)
    (let
      in do sourceR2S1RP <-
              r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs maxScale $
              inputR2Z2
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
  -> Double
  -> R2Z2T0S0Array
  -> PowerMethodNormalizationOption
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
powerMethodR2Z2T0S0Bias plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs maxScale filter normMethod numIteration writeFlag idStr threshold bias eigenVecSource = do
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
           maxScale
           filterF
           normMethod
           n
           writeFlag
           ("Source" L.++ idStr)
           bias
           input)
      (computeS $ delay eigenVecSource)
      [0 .. numIteration]
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

{-# INLINE eigenVectorR2Z2Reversal #-}
eigenVectorR2Z2Reversal
  --    (R.Source s1 (Complex Double), R.Source s2 (Complex Double))
  -- =>
 ::
     DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> Double
  -> R2Z2T0S0Array
  -> PowerMethodNormalizationOption
  -> Int
  -> Bool
  -> String
  -> Double
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> (R.Array U DIM4 (Complex Double), R.Array U DIM6 (Complex Double))
  -> IO (R2T0S0Array, R2Z2T0S0Array)
eigenVectorR2Z2Reversal plan folderPath numOrientation thetaFreqs numScale scaleFreqs maxScale filterF normMethod n writeFlag name reversalFactor approximatedEigenValue bias (inputR2Z2, _inputR2Z2T0S0)
      -- (Z :. numTheta0Freq :. numScale0Freq :. cols :. rows :. numThetaFreq :. numScaleFreq) =
      --   extent inputR2Z2T0S0
      -- (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
      --   extent inputR2Z2T0S0
 = do
  let (Z :. numThetaFreq :. numScaleFreq :. cols :. rows) = extent inputR2Z2
  -- inputR2Z2 <-
  --   (R.sumP .
  --    rotateR2Z2T0S0ArrayDeconv -- . computeSinkFromSourceR2Z2T0S0 thetaFreqs thetaFreqs
  --    $
  --    inputR2Z2T0S0) >>=
  --   R.sumP
  let biasedInputR2Z2 = R.zipWith (*) bias inputR2Z2
  -- s <- R.sumAllP . rotate4D . rotate4D . R.map magnitude $ biasedInputR2Z2
  normalizedBiasedInputR2Z2 <-
    powerMethodNormalization normMethod biasedInputR2Z2
  -- let normalizedBiasedInputR2Z2 = computeUnboxedS . R.map (/ (s :+ 0)) $ inputR2Z2
  let initialDist =
        R.traverse2
          normalizedBiasedInputR2Z2
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          const $ \f1 f2 idx@(Z :. k :. l :. i :. j) ->
          f1 idx + (reversalFactor :+ 0) * f1 idx * exp (0 :+ f2 (Z :. k) * pi)
  -- let initialDist =
  --       R.traverse2
  --         normalizedBiasedInputR2Z2
  --         (R.sumS . R.sumS . R.map magnitude . rotate4D . rotate4D $
  --          normalizedBiasedInputR2Z2)
  --         const $ \f1 f2 idx@(Z :. k :. l :. i :. j) ->
  --         if k == div numThetaFreq 2 && l == div numScaleFreq 2
  --           then f1 idx + (reversalFactor * f2 (Z :. i :. j) :+ 0)
  --           else f1 idx
  sourceArr <- convolveR2T0S0P plan filterF . computeUnboxedS $ initialDist
  sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceArr) >>= R.sumP -- >>=
    -- computeP .
    -- R.zipWith (\x y -> y - (approximatedEigenValue :+ 0) * x ) initialDist
  let a =
        sumAllS $
        R.zipWith
          (\x y -> conjugate x * y)
          initialDist
          (R.zipWith (*) bias sourceR2Z2)
      b = sumAllS $ R.zipWith (\x y -> conjugate x * y) initialDist initialDist
  printCurrentTime $
    printf "iteration %d: %.5f" (n + 1) ((magnitude a) / (magnitude b))
  when
    (n == 0 || writeFlag)
    (let
      in do sourceR2S1RP <-
              r2z2Tor2s1rpP
                numOrientation
                thetaFreqs
                numScale
                scaleFreqs
                maxScale $
              sourceR2Z2
            sourceField <-
              (R.sumP . rotate4D . rotate4D $ sourceR2S1RP) >>=
              fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
              R.sumP
            plotImageRepaComplex
              (folderPath </> name L.++ "_" L.++ show n L.++ ".png") .
              ImageRepa 8 $
              sourceField
            let biasR2Z2 = R.zipWith (*) bias sourceR2Z2
            biasR2S1RP <-
              r2z2Tor2s1rpP
                numOrientation
                thetaFreqs
                numScale
                scaleFreqs
                maxScale $
              biasR2Z2
            biasField <-
              (R.sumP . rotate4D . rotate4D . R.map magnitude $ biasR2S1RP) >>=
              fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
              R.sumP
            plotImageRepa (folderPath </> "Bias_" L.++ show n L.++ ".png") .
              ImageRepa 8 $
              biasField
            -- let nonzeroVec = VU.filter (/= 0) . toUnboxed $ biasField
            --     minV = VU.minimum nonzeroVec
            --     maxV = VU.maximum nonzeroVec
            -- plotImageRepa
            --   (folderPath </> "InversedBias_" L.++ show n L.++ ".png") .
            --   ImageRepa 8 .
            --   computeS .
            --   R.map
            --     (\x ->
            --        if x == 0
            --          then 0
            --          else (maxV - x) + minV) $
            --   biasField
     )
  return (sourceR2Z2, sourceArr)

powerMethodR2Z2T0S0Reversal ::
     (R.Source s2 (Complex Double))
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
  -> Double
  -> R2Z2T0S0Array
  -> PowerMethodNormalizationOption
  -> Int
  -> Bool
  -> String
  -> Double
  -> Double
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s2 DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0Reversal plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs maxScale filter normMethod numIteration writeFlag idStr threshold reversalFactor approximatedEigenValue bias eigenVecSource = do
  filterF <-
    dftR2Z2T0S0 plan .
    computeS .
    makeFilter2D .
    -- computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs .
    -- R.traverse filter id $ \f idx@(Z :. a :. b :. c :. d :. i :. j) ->
    --   if a == div (L.length thetaFreqs) 2 &&
    --      b == div (L.length scaleFreqs) 2 &&
    --      c == div (L.length theta0Freqs) 2 && d == div (L.length scale0Freqs) 2
    --     then 1 +
    --          conjugate
    --            (f (Z :. (L.length thetaFreqs - 1 - a) :.
    --                (L.length scaleFreqs - 1 - b) :.
    --                (L.length theta0Freqs - 1 - c) :.
    --                (L.length scale0Freqs - 1 - d) :.
    --                i :.
    --                j))
    --     else conjugate
    --            (f (Z :. (L.length thetaFreqs - 1 - a) :.
    --                (L.length scaleFreqs - 1 - b) :.
    --                (L.length theta0Freqs - 1 - c) :.
    --                (L.length scale0Freqs - 1 - d) :.
    --                i :.
    --                j))
    R.traverse (computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs filter) id $ \f (Z :. tf' :. sf' :. t0f :. s0f :. i :. j) ->
      let idx = (Z :. tf' :. sf' :. t0f :. s0f :. i :. j)
       in if i == center cols && j == center rows && tf' == t0f
             -- tf' == div (L.length thetaFreqs) 2 &&
             -- sf' == div (L.length scaleFreqs) 2 &&
             -- t0f == div (L.length thetaFreqs) 2 &&
             -- s0f == div (L.length scaleFreqs) 2
            then f idx - (approximatedEigenValue :+ 0)
            else f idx
  -- filterF <-
  --   dftR2Z2T0S0Deconv plan .
  --   computeS . makeFilter4D . rotateR2Z2T0S0ArrayDeconv . R.traverse filter id $ \f (Z :. tf' :. sf' :. t0f :. s0f :. i :. j) ->
  --     let -- idx =
  --         --   (Z :. ((L.length thetaFreqs) - 1 - tf') :.
  --         --    ((L.length scaleFreqs) - 1 - sf') :.
  --         --    t0f :.
  --         --    s0f :.
  --         --    i :.
  --         --    j)
  --         idx =
  --           (Z :. tf' :.
  --            sf' :.
  --            t0f :.
  --            s0f :.
  --            i :.
  --            j)
  --      in if i == center cols && j == center rows
  --                                                                                                                       -- t0f == div (L.length thetaFreqs) 2 && s0f == div (L.length scaleFreqs) 2
  --           then f idx - (approximatedEigenValue :+ 0)
  --           else f idx
  initR2Z2 <-
    (R.sumP .
     rotateR2Z2T0S0Array -- . computeSinkFromSourceR2Z2T0S0 thetaFreqs thetaFreqs
      $
     eigenVecSource) >>=
    R.sumP
  (sourceR2Z2, sourceR2Z2T0S0) <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2Reversal
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           maxScale
           filterF
           normMethod
           n
           writeFlag
           ("Source" L.++ idStr)
           reversalFactor
           approximatedEigenValue
           bias
           input)
      (initR2Z2, computeS . delay $ eigenVecSource)
      [0 .. numIteration]
  let sinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs sourceR2Z2T0S0
  -- sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0ArrayDeconv $ sourceR2Z2T0S0) >>= R.sumP
  sinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sinkR2Z2T0S0) >>= R.sumP
  sinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     sinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s.png" idStr) . ImageRepa 8 $
    sinkField
  -- R.zipWith (*) bias <$>
  _ <-
    delay <$>
    completionFieldR2Z2
      plan
      folderPath
      idStr
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      (computeS . R.zipWith (\x y -> (1-x) * y) bias $ sourceR2Z2)
      (computeS . R.zipWith (\x y -> (1-x) * y) bias $ sinkR2Z2)
  return . delay $ sourceR2Z2

{-# INLINE eigenVectorR2Z2BiasReversal #-}
eigenVectorR2Z2BiasReversal ::
     (R.Source s (Complex Double))
  => DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> Double
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO R2Z2T0S0Array
eigenVectorR2Z2BiasReversal plan folderPath numOrientation thetaFreqs numScale scaleFreqs maxScale filterF n writeFlag name reversalFactor bias inputR2Z2T0S0 = do
  let (Z :. numThetaFreq :. numScaleFreq :. numTheta0Freq :. numScale0Freq :. cols :. rows) =
        extent inputR2Z2T0S0
  inputR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ inputR2Z2T0S0) >>= R.sumP
  sourceDist <- convolveR2Z2 plan bias inputR2Z2
  s <- R.foldAllP max 0 . R.map magnitude $ sourceDist
  let initialDist' = R.map (/ (s :+ 0)) sourceDist
      initialDist =
        R.zipWith (\x y -> x + (reversalFactor :+ 0) * y) initialDist' .
        R.traverse2
          initialDist'
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          const $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
          f1 idx * exp (0 :+ f2 (Z :. k) * pi)
  sourceArr <- convolveR2T0S0P plan filterF initialDist
  printCurrentTime $ printf "iteration %d" n
  when
    (n == 1 ||
     (writeFlag -- && odd n
      ))
    (let
      in do sourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ sourceArr) >>= R.sumP
            sourceR2S1RP <-
              r2z2Tor2s1rpP
                numOrientation
                thetaFreqs
                numScale
                scaleFreqs
                maxScale $
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

powerMethodR2Z2T0S0BiasReversal ::
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
  -> Double
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0BiasReversal plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs maxScale filter numIteration writeFlag idStr reversalFactor bias eigenVecSource = do
  filterF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filter
  sourceR2Z2T0S0 <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2BiasReversal
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           maxScale
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
  -- R.zipWith (*) bias <$>
  delay <$>
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



{-# INLINE eigenVectorR2Z2EndModal #-}
eigenVectorR2Z2EndModal ::
     DFTPlan
  -> FilePath
  -> Int
  -> [Double]
  -> Int
  -> [Double]
  -> Double
  -> R2Z2T0S0Array
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> (R.Array U DIM4 (Complex Double), R.Array U DIM4 (Complex Double))
  -> IO (R2T0S0Array, R2T0S0Array)
eigenVectorR2Z2EndModal plan folderPath numOrientation thetaFreqs numScale scaleFreqs maxScale filterEndF filterModalF n writeFlag name reversalFactor bias (endR2Z2, modalR2Z2) = do
  let endBias =
        computeS $
        R.zipWith
          (+)
          (rotateBiasR2Z2T0S0 90 thetaFreqs modalR2Z2)
          (rotateBiasR2Z2T0S0 (-90) thetaFreqs modalR2Z2)
  endEigenVec <-
    M.foldM
      (\input _
         -- endInitialDist' <- convolveR2Z2 plan endBias input
        -> do
         let endInitialDist =
               R.traverse2
                 -- endInitialDist'
                 input
                 (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
                 const $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
                 f1 idx +
                 (reversalFactor :+ 0) * f1 idx * exp (0 :+ f2 (Z :. k) * pi)
         endSourceArr <- convolveR2T0S0P plan filterEndF endInitialDist
         endSourceR2Z2 <-
           (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
         biasedEndSourceR2Z2 -- convolveR2Z2 plan bias endSourceR2Z2
            <-
           convolveR2Z2 plan endBias endSourceR2Z2 >>= convolveR2Z2 plan bias
         sEnd <- R.foldAllP max 0 . R.map magnitude $ biasedEndSourceR2Z2
         return . computeS $ R.map (/ (sEnd :+ 0)) biasedEndSourceR2Z2)
      endR2Z2
      [1 .. 2]
  -- endInitialDist' <- convolveR2Z2 plan endBias endEigenVec
  let endInitialDist =
        R.traverse2
          -- endInitialDist'
          endEigenVec
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          const $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
          f1 idx + (reversalFactor :+ 0) * f1 idx * exp (0 :+ f2 (Z :. k) * pi)
  endSourceArr <- convolveR2T0S0P plan filterEndF endInitialDist
  endSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
  let modalBias =
        computeS $
        R.zipWith
          (+)
          (rotateBiasR2Z2T0S0 90 thetaFreqs endEigenVec)
          (rotateBiasR2Z2T0S0 (-90) thetaFreqs endEigenVec)
  modalEigenVec <-
    M.foldM
      (\input _
         -- modalInitialDist <- convolveR2Z2 plan modalBias input
        -> do
         modalSourceArr <- convolveR2T0S0P plan filterModalF input -- modalInitialDist
         modalSourceR2Z2 <-
           (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
         biasedModalSourceR2Z2 -- convolveR2Z2 plan bias modalSourceR2Z2
            <-
           convolveR2Z2 plan modalBias modalSourceR2Z2 >>=
           convolveR2Z2 plan bias
         sModal <- R.foldAllP max 0 . R.map magnitude $ biasedModalSourceR2Z2
         return . computeS $ R.map (/ (sModal :+ 0)) biasedModalSourceR2Z2)
      modalR2Z2
      [1 .. 2]
  -- modalInitialDist <- convolveR2Z2 plan modalBias modalEigenVec
  modalSourceArr <- convolveR2T0S0P plan filterModalF modalEigenVec --modalInitialDist
  modalSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
  printCurrentTime $ printf "iteration %d" n
  when
    writeFlag
    (do endSourceR2S1RP <-
          r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs maxScale $
          endSourceR2Z2
        endSourceField <-
          (R.sumP . rotate4D . rotate4D . R.map magnitude $ endSourceR2S1RP) >>=
          fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
        plotImageRepa
          (folderPath </> name L.++ "_Endpoint_" L.++ show n L.++ ".png") .
          ImageRepa 8 $
          endSourceField
        modalSourceR2Z2 <-
          (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
        modalSourceR2S1RP <-
          r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs maxScale $
          modalSourceR2Z2
        modalSourceField <-
          (R.sumP . rotate4D . rotate4D . R.map magnitude $ modalSourceR2S1RP) >>=
          fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
        plotImageRepa
          (folderPath </> name L.++ "_Modal_" L.++ show n L.++ ".png") .
          ImageRepa 8 $
          modalSourceField)
  return (endEigenVec, modalEigenVec)


powerMethodR2Z2T0S0EndModal ::
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
  -> Double
  -> R2Z2T0S0Array
  -> R2Z2T0S0Array
  -> Int
  -> Bool
  -> String
  -> Double
  -> R.Array U DIM4 (Complex Double)
  -> R.Array s DIM6 (Complex Double)
  -> IO (R.Array D DIM4 (Complex Double), R.Array D DIM4 (Complex Double))
powerMethodR2Z2T0S0EndModal plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs maxScale filterEnd filterModal numIteration writeFlag idStr reversalFactor bias eigenVecSource = do
  filterEndF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filterEnd
  filterModalF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filterModal
  endR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ eigenVecSource) >>= R.sumP
  endInitialDist' <- convolveR2Z2 plan bias endR2Z2
  let endInitialDist =
        R.traverse2
          endInitialDist'
          (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
          const $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
          f1 idx + (reversalFactor :+ 0) * f1 idx * exp (0 :+ f2 (Z :. k) * pi)
  endSourceArr <- convolveR2T0S0P plan filterEndF endInitialDist
  endSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
  biasedEndSourceR2Z2 <- convolveR2Z2 plan bias endSourceR2Z2
  sEnd <- R.foldAllP max 0 . R.map magnitude $ biasedEndSourceR2Z2
  let endEigenVec = computeS $ R.map (/ (sEnd :+ 0)) biasedEndSourceR2Z2
      modalBias =
        computeS $
        R.zipWith
          (+)
          (rotateBiasR2Z2T0S0 90 thetaFreqs endEigenVec)
          (rotateBiasR2Z2T0S0 (-90) thetaFreqs endEigenVec)
  modalInitialDist <- convolveR2Z2 plan modalBias endR2Z2
  modalSourceArr <- convolveR2T0S0P plan filterModalF modalInitialDist
  modalSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
  biasedModalSourceR2Z2 <-
    convolveR2Z2 plan modalBias modalSourceR2Z2 >>= convolveR2Z2 plan bias
  sModal <- R.foldAllP max 0 . R.map magnitude $ biasedModalSourceR2Z2
  let modalEigenVec = computeS $ R.map (/ (sModal :+ 0)) biasedModalSourceR2Z2
  -- plot 
  endSourceR2S1RP <-
    r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs maxScale $ endSourceR2Z2
  endSourceField <-
    (R.sumP . rotate4D . rotate4D $ endSourceR2S1RP) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex
    (folderPath </> "Source" L.++ idStr L.++ "_Endpoint_0.png") .
    ImageRepa 8 $
    endSourceField
  modalSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
  modalSourceR2S1RP <-
    r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs maxScale $
    modalSourceR2Z2
  modalSourceField <-
    (R.sumP . rotate4D . rotate4D $ modalSourceR2S1RP) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> "Source" L.++ idStr L.++ "_Modal_0.png") .
    ImageRepa 8 $
    modalSourceField
  --
  (endEigenVecR2Z2, modalEigenVecR2Z2) <-
    M.foldM
      (\input n ->
         eigenVectorR2Z2EndModal
           plan
           folderPath
           numOrientation
           thetaFreqs
           numScale
           scaleFreqs
           maxScale
           filterEndF
           filterModalF
           n
           writeFlag
           ("Source" L.++ idStr)
           reversalFactor
           bias
           input)
      (endEigenVec, modalEigenVec)
      [1 .. numIteration]
  modalSourceR2Z2T0S0 <- convolveR2T0S0P plan filterModalF modalEigenVecR2Z2
  endSourceR2Z2T0S0 <- convolveR2T0S0P plan filterEndF endEigenVecR2Z2
  let modalSinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs modalSourceR2Z2T0S0
  modalSourceR2Z2 <-
    (R.sumP . rotateR2Z2T0S0Array $ modalSourceR2Z2T0S0) >>= R.sumP
  modalSinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalSinkR2Z2T0S0) >>= R.sumP
  modalSinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     modalSinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s_Modal.png" idStr) .
    ImageRepa 8 $
    modalSinkField
  let endSinkR2Z2T0S0 =
        computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs endSourceR2Z2T0S0
  endSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceR2Z2T0S0) >>= R.sumP
  endSinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSinkR2Z2T0S0) >>= R.sumP
  endSinkField <-
    (R.sumP .
     rotate4D .
     rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
     endSinkR2Z2) >>=
    fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
  plotImageRepaComplex (folderPath </> printf "Sink%s.png_EndPoint" idStr) .
    ImageRepa 8 $
    endSinkField
  completionModal <-
    delay <$>
    completionFieldR2Z2
      plan
      folderPath
      (idStr L.++ "_Modal")
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      modalSourceR2Z2
      modalSinkR2Z2
  completionEnd <-
    delay <$>
    completionFieldR2Z2
      plan
      folderPath
      (idStr L.++ "_EndPoint")
      numOrientation
      thetaFreqs
      numScale
      scaleFreqs
      endSourceR2Z2
      endSinkR2Z2
  return (completionEnd, completionModal)



-- {-# INLINE eigenVectorR2Z2EndModal #-}
-- eigenVectorR2Z2EndModal ::
--      (R.Source s (Complex Double))
--   => DFTPlan
--   -> FilePath
--   -> Int
--   -> [Double]
--   -> Int
--   -> [Double]
--   -> R2Z2T0S0Array
--   -> R2Z2T0S0Array
--   -> Int
--   -> Bool
--   -> String
--   -> Double
--   -> R.Array U DIM4 (Complex Double)
--   -> (R.Array s DIM6 (Complex Double), R.Array s DIM6 (Complex Double))
--   -> IO (R2Z2T0S0Array, R2Z2T0S0Array)
-- eigenVectorR2Z2EndModal plan folderPath numOrientation thetaFreqs numScale scaleFreqs filterEndF filterModalF n writeFlag name reversalFactor bias (endR2Z2T0S0, modalR2Z2T0S0) = do
--   endR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endR2Z2T0S0) >>= R.sumP
--   modalR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalR2Z2T0S0) >>= R.sumP
--   let endBias =
--         computeS $
--         R.zipWith
--           (+)
--           (rotateBiasR2Z2T0S0 90 thetaFreqs modalR2Z2)
--           (rotateBiasR2Z2T0S0 (-90) thetaFreqs modalR2Z2)
--   endSourceDist' <- convolveR2Z2 plan endBias endR2Z2
--   endSourceDist <- convolveR2Z2 plan bias endSourceDist'
--   sEnd <- R.foldAllP max 0 . R.map magnitude $ endSourceDist
--   let endInitialDist' = R.map (/ (sEnd :+ 0)) endSourceDist
--       endInitialDist =
--         R.traverse2
--           endInitialDist'
--           (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
--           const $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
--           f1 idx + (reversalFactor :+ 0) * f1 idx * exp (0 :+ f2 (Z :. k) * pi)
--   endSourceArr <- convolveR2T0S0P plan filterEndF endInitialDist
--   newEndR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
--   let modalBias =
--         computeS $
--         R.zipWith
--           (+)
--           (rotateBiasR2Z2T0S0 90 thetaFreqs newEndR2Z2)
--           (rotateBiasR2Z2T0S0 (-90) thetaFreqs newEndR2Z2)
--   modalSourceDist' <- convolveR2Z2 plan modalBias modalR2Z2
--   modalSourceDist <- convolveR2Z2 plan bias modalSourceDist'
--   sModal <- R.foldAllP max 0 . R.map magnitude $ modalSourceDist
--   let modalInitialDist = R.map (/ (sModal :+ 0)) modalSourceDist
--   modalSourceArr <- convolveR2T0S0P plan filterModalF modalInitialDist
--   printCurrentTime $ printf "iteration %d" n
--   when
--     (writeFlag -- && odd n
--      )
--     (let
--       in do endSourceR2Z2 <-
--               (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
--             endSourceR2S1RP <-
--               r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
--               endSourceR2Z2
--             endSourceField <-
--               (R.sumP . rotate4D . rotate4D $ endSourceR2S1RP) >>=
--               fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
--               R.sumP
--             plotImageRepaComplex
--               (folderPath </> name L.++ "_Endpoint_" L.++ show n L.++ ".png") .
--               ImageRepa 8 $
--               endSourceField
--             modalSourceR2Z2 <-
--               (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
--             modalSourceR2S1RP <-
--               r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
--               modalSourceR2Z2
--             modalSourceField <-
--               (R.sumP . rotate4D . rotate4D $ modalSourceR2S1RP) >>=
--               fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) .
--               R.sumP
--             plotImageRepaComplex
--               (folderPath </> name L.++ "_Modal_" L.++ show n L.++ ".png") .
--               ImageRepa 8 $
--               modalSourceField)
--   return (endSourceArr, modalSourceArr)


-- powerMethodR2Z2T0S0EndModal ::
--      (R.Source s (Complex Double))
--   => DFTPlan
--   -> FilePath
--   -> Int
--   -> Int
--   -> Int
--   -> [Double]
--   -> [Double]
--   -> Int
--   -> [Double]
--   -> [Double]
--   -> R2Z2T0S0Array
--   -> R2Z2T0S0Array
--   -> Int
--   -> Bool
--   -> String
--   -> Double
--   -> R.Array U DIM4 (Complex Double)
--   -> R.Array s DIM6 (Complex Double)
--   -> IO (R.Array D DIM4 (Complex Double), R.Array D DIM4 (Complex Double))
-- powerMethodR2Z2T0S0EndModal plan folderPath cols rows numOrientation thetaFreqs theta0Freqs numScale scaleFreqs scale0Freqs filterEnd filterModal numIteration writeFlag idStr reversalFactor bias eigenVecSource = do
--   filterEndF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filterEnd
--   filterModalF <- dftR2Z2T0S0 plan . computeS . makeFilter2D $ filterModal
--   endR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ eigenVecSource) >>= R.sumP
--   endSourceDist <- convolveR2Z2 plan bias endR2Z2
--   sEnd <- R.foldAllP max 0 . R.map magnitude $ endSourceDist
--   let endInitialDist' = R.map (/ (sEnd :+ 0)) endSourceDist
--       endInitialDist =
--         R.traverse2
--           endInitialDist'
--           (fromListUnboxed (Z :. L.length thetaFreqs) thetaFreqs)
--           const $ \f1 f2 idx@(Z :. k :. _ :. i :. j) ->
--           f1 idx + (reversalFactor :+ 0) * f1 idx * exp (0 :+ f2 (Z :. k) * pi)
--   endSourceArr <- convolveR2T0S0P plan filterEndF endInitialDist
--   newEndR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
--   let modalBias =
--         computeS $
--         R.zipWith
--           (+)
--           (rotateBiasR2Z2T0S0 90 thetaFreqs newEndR2Z2)
--           (rotateBiasR2Z2T0S0 (-90) thetaFreqs newEndR2Z2)
--   modalSourceDist' <- convolveR2Z2 plan modalBias endR2Z2
--   modalSourceDist <- convolveR2Z2 plan bias modalSourceDist'
--   sModal <- R.foldAllP max 0 . R.map magnitude $ modalSourceDist
--   let modalInitialDist = R.map (/ (sModal :+ 0)) modalSourceDist
--   modalSourceArr <- convolveR2T0S0P plan filterModalF modalInitialDist
--   endSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceArr) >>= R.sumP
--   endSourceR2S1RP <-
--     r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $ endSourceR2Z2
--   endSourceField <-
--     (R.sumP . rotate4D . rotate4D $ endSourceR2S1RP) >>=
--     fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
--   plotImageRepaComplex
--     (folderPath </> "Source" L.++ idStr L.++ "_Endpoint_0.png") .
--     ImageRepa 8 $
--     endSourceField
--   modalSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalSourceArr) >>= R.sumP
--   modalSourceR2S1RP <-
--     r2z2Tor2s1rpP numOrientation thetaFreqs numScale scaleFreqs $
--     modalSourceR2Z2
--   modalSourceField <-
--     (R.sumP . rotate4D . rotate4D $ modalSourceR2S1RP) >>=
--     fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
--   plotImageRepaComplex (folderPath </> "Source" L.++ idStr L.++ "_Modal_0.png") .
--     ImageRepa 8 $
--     modalSourceField
--   (endSourceR2Z2T0S0, modalSourceR2Z2T0S0) <-
--     M.foldM
--       (\input n ->
--          eigenVectorR2Z2EndModal
--            plan
--            folderPath
--            numOrientation
--            thetaFreqs
--            numScale
--            scaleFreqs
--            filterEndF
--            filterModalF
--            n
--            writeFlag
--            ("Source" L.++ idStr)
--            reversalFactor
--            bias
--            input)
--       (endSourceArr, modalSourceArr)
--       [1 .. numIteration]
--   let modalSinkR2Z2T0S0 =
--         computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs modalSourceR2Z2T0S0
--   modalSourceR2Z2 <-
--     (R.sumP . rotateR2Z2T0S0Array $ modalSourceR2Z2T0S0) >>= R.sumP
--   modalSinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ modalSinkR2Z2T0S0) >>= R.sumP
--   modalSinkField <-
--     (R.sumP .
--      rotate4D .
--      rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
--      modalSinkR2Z2) >>=
--     fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
--   plotImageRepaComplex (folderPath </> printf "Sink%s_Modal.png" idStr) .
--     ImageRepa 8 $
--     modalSinkField
--   let endSinkR2Z2T0S0 =
--         computeSinkFromSourceR2Z2T0S0 thetaFreqs theta0Freqs endSourceR2Z2T0S0
--   endSourceR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSourceR2Z2T0S0) >>= R.sumP
--   endSinkR2Z2 <- (R.sumP . rotateR2Z2T0S0Array $ endSinkR2Z2T0S0) >>= R.sumP
--   endSinkField <-
--     (R.sumP .
--      rotate4D .
--      rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
--      endSinkR2Z2) >>=
--     fmap (computeS . R.extend (Z :. (1 :: Int) :. All :. All)) . R.sumP
--   plotImageRepaComplex (folderPath </> printf "Sink%s.png_EndPoint" idStr) .
--     ImageRepa 8 $
--     endSinkField
--   completionModal <-
--     delay <$>
--     completionFieldR2Z2
--       plan
--       folderPath
--       (idStr L.++ "_Modal")
--       numOrientation
--       thetaFreqs
--       numScale
--       scaleFreqs
--       modalSourceR2Z2
--       modalSinkR2Z2
--   completionEnd <-
--     delay <$>
--     completionFieldR2Z2
--       plan
--       folderPath
--       (idStr L.++ "_EndPoint")
--       numOrientation
--       thetaFreqs
--       numScale
--       scaleFreqs
--       endSourceR2Z2
--       endSinkR2Z2
--   return (completionEnd, completionModal)

