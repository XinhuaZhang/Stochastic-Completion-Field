{-# LANGUAGE BangPatterns #-}
module LogFourierSeries where

import           Data.Array.Repa     as R
import           Data.Complex
import           Data.List           as L
import           Data.Vector.Unboxed as VU
import           Image.IO
import           System.Environment
import           System.FilePath
import           Utils.Parallel hiding ((.|))
import Control.DeepSeq
import Data.Conduit
import Data.Conduit.List as CL
import Control.Monad.Trans.Resource
import Control.Monad
import           FokkerPlanck.MonteCarlo
import Utils.Array
import            STC.Utils



{-# INLINE gaussian #-}
gaussian :: Double -> Double -> Double -> Double
gaussian !mu !sigma !x =
  (exp $ (x - mu) ^ 2 / (-2 * sigma ^ 2)) / (sqrt $ 2 * pi) / sigma
  
sink ::
     ParallelParams -> (VU.Vector (Complex Double))
  -> ConduitT (VU.Vector (Complex Double)) Void (ResourceT IO) (VU.Vector Double)
sink !parallelParams !vec = do
  xs <- CL.take . batchSize $ parallelParams
  if (L.null xs)
    then return . VU.map magnitude $ vec
    else sink parallelParams (L.foldl' (VU.zipWith (+)) vec xs)
  -- case x of
  --   Nothing -> return . VU.map magnitude $ vec
  --   Just y -> sink (VU.zipWith (+) vec y)

main = do
  args@(numPointStr:maxThetaFreqStr:maxRFreqStr:deltaStr:maxRStr:stdStr:scaleFactorStr:batchSizeStr:numThreadStr:_) <-
    getArgs
  print args
  let maxThetaFreq = read maxThetaFreqStr :: Double
      maxRFreq = read maxRFreqStr :: Double
      delta = read deltaStr :: Double
      thetaFreqs = [-maxThetaFreq .. maxThetaFreq]
      rFreqs = [-maxRFreq .. maxRFreq]
      phiFreqs = [0] -- thetaFreqs
      muR = 16
      sigmaR = 8
      muTheta = 0
      sigmaTheta = 0.5
      muPhi = pi
      sigmaPhi = 0.5
      deltaPhi = 2 * pi / fromIntegral numOrientation
      alpha = 1
      scaleFactor = read scaleFactorStr :: Double
      numPoint = read numPointStr :: Int
      numOrientation = 72 :: Int
      maxR = read maxRStr :: Double
      std = read stdStr :: Double
      batchSize = read batchSizeStr :: Int
      numThread = read numThreadStr :: Int
      center = div numPoint 2
      folderPath = "output/test/LogFourierSeries"
      gaussianFunc theta r phi =
        (gaussian muTheta sigmaTheta theta) * (gaussian muR sigmaR r) *
        (gaussian muPhi sigmaPhi phi)
      centerArr =
        fromFunction (Z :. (1 :: Int) :. numPoint :. numPoint) $ \(Z :. _ :. i :. j) ->
          if i == center && j == center
            then 1
            else 0
  plotImageRepa (folderPath </> "Center.png") . ImageRepa 8 . computeS $
    centerArr
  -- originalArray <-
  --   computeUnboxedP . fromFunction (Z :. numPoint :. numPoint :. numOrientation) $ \(Z :. i :. j :. k) ->
  --     let !x = (fromIntegral $ i)
  --         !y = (fromIntegral $ j - center)
  --         !theta = atan2 y x
  --         !r = sqrt $ x ^ 2 + y ^ 2
  --         !phi = fromIntegral k * deltaPhi
  --     in gaussianFunc theta r phi
  arrG <-
    solveMonteCarloR2S1
      numThread
      1000000
      1000000
      numPoint
      numPoint
      numOrientation
      0.1
      100
      64
      1
      ""
  print . extent $ arrG
  let originalArray =
        extend (Z :. All :. All :. (1 :: Int)) . sumS $ rotate3D arrG
  plotImageRepa (folderPath </> "Greens3DNorm_1.png") .
    ImageRepa 8 .
    computeS .
    reduceContrast 10 . extend (Z :. (1 :: Int) :. All :. All) . R.sumS $
    originalArray
  plotImageRepa (folderPath </> "Greens3D_1.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) . R.sumS $
    originalArray
  recon <-
    runConduitRes $
    CL.sourceList
      [ (rFreq, thetaFreq, phiFreq)
      | rFreq <- rFreqs
      , thetaFreq <- thetaFreqs
      , phiFreq <- phiFreqs
      ] .|
    parConduit
      (ParallelParams numThread batchSize)
      (\(rFreq, thetaFreq, phiFreq) ->
         let !c =
               sumAllS .
               R.zipWith (\a b -> (a :+ 0) * b) originalArray .
               fromFunction (Z :. numPoint :. numPoint :. (1 :: Int)) $ \(Z :. i :. j :. k) ->
                 let !x = delta * (fromIntegral $ i - center)
                     !y = delta * (fromIntegral $ j - center)
                     !theta = atan2 y x
                     !r = (sqrt $ x ^ 2 + y ^ 2)
                     !phi = fromIntegral k * deltaPhi
                 in if (x ^ 2 + y ^ 2) == 0  -- || r >= maxR -- || pi * r <  (abs thetaFreq)
                      then 0
                           -- exp $ (log r) * (alpha - 1) :+ ((-thetaFreq) * theta - rFreq * (log r))
                      else exp $ (log r) * (-1 + alpha) :+ ((-thetaFreq) * theta - rFreq * (log r) -- * pi / (log maxR)
                                                              )
                           -- cis ((-thetaFreq) * theta - rFreq * (log r) * pi / (log maxR) )
                           -- (x :+ y) ** ((-thetaFreq) :+ 0) *
                           -- ((x ^ 2 + y ^ 2) :+ 0) **
                           -- (((thetaFreq - 2 + alpha) :+ (-rFreq)) / 2)
                           -- (x :+ y) ** ((-thetaFreq) :+ 0) *
                           -- ((x ^ 2 + y ^ 2) :+ 0) **
                           -- (((thetaFreq) :+ (-pi * rFreq / (log maxR))) / 2)
         in toUnboxed .
            computeS .
            -- R.map (* c) .
            R.map
              (* (c *
                  ((exp $ (-thetaFreq ^ 2 - rFreq ^ 2) / (2 * std ^ 2)) :+ 0))) .
            fromFunction (Z :. numPoint :. numPoint :. (1 :: Int)) $ \(Z :. i :. j :. k) ->
              let !x = (fromIntegral $ i - center)
                  !y = (fromIntegral $ j - center)
                  !theta = atan2 y x
                  !r = delta * (sqrt $ x ^ 2 + y ^ 2)
                  !phi = fromIntegral k * deltaPhi
              in if r <= 0 -- || pi * r < (abs thetaFreq)
                   then 0
                        
                   else exp $ (log r) * (-alpha) :+ (thetaFreq * theta + rFreq * (log r) -- * pi / (log maxR)
                                                  )
                        -- cis (thetaFreq * theta + rFreq * (log r) * pi / (log maxR))
                        -- exp $ (log r) * (-alpha) :+ (thetaFreq * theta + rFreq * (log r))
                        -- (x :+ y) ** (thetaFreq :+ 0) *
                        -- ((x ^ 2 + y ^ 2) :+ 0) **
                        -- (((-thetaFreq - alpha) :+ rFreq) / 2)
                        -- (x :+ y) ** (thetaFreq :+ 0) *
                        -- ((x ^ 2 + y ^ 2) :+ 0) **
                        -- (((-thetaFreq) :+ (pi * rFreq / (log maxR))) / 2)
       ) .|
    sink
      (ParallelParams numThread batchSize)
      (VU.replicate (numPoint ^ 2 * 1) 0)
  -- let
      -- !recon =
      --   VU.map magnitude .
      --   L.foldl1' (VU.zipWith (+)) .
      --   parMap
      --     rdeepseq
      --     (\(rFreq, thetaFreq, phiFreq) ->
      --        let !c =
      --              sumAllS .
      --              fromFunction (Z :. numPoint :. numPoint :. numPoint) $ \(Z :. i :. j :. k) ->
      --                let !x = delta * (fromIntegral $ i - center)
      --                    !y = delta * (fromIntegral $ j - center)
      --                    !theta = atan2 y x
      --                    !r = sqrt $ x ^ 2 + y ^ 2
      --                    !phi = fromIntegral k * deltaPhi
      --                in if r < 1 / maxR || r >= maxR
      --                     then 0
      --                     else (gaussianFunc theta r phi :+ 0) *
      --                          (cis $
      --                           (-1) * thetaFreq * theta -
      --                           phiFreq * (phi - theta) -
      --                           pi * rFreq * (r) / (maxR))
      --        in toUnboxed .
      --           computeS .
      --           R.map (* c) .
      --           fromFunction (Z :. numPoint :. numPoint :. numPoint) $ \(Z :. i :. j :. k) ->
      --             let !x = fromIntegral $ i - center
      --                 !y = fromIntegral $ j - center
      --                 !theta = atan2 y x
      --                 !r = sqrt $ x ^ 2 + y ^ 2
      --                 !phi = fromIntegral k * deltaPhi
      --             in if r == 0
      --                  then 0
      --                  else cis $
      --                       thetaFreq * theta + phiFreq * (phi - theta) +
      --                       pi * rFreq * (r) / (maxR)) $
      --   [ (rFreq, thetaFreq, phiFreq)
      --   | rFreq <- rFreqs
      --   , thetaFreq <- thetaFreqs
      --   , phiFreq <- phiFreqs
      --   ]
  plotImageRepa (folderPath </> "Greens3DRecon_2.png") .
    ImageRepa 8 .
    computeS .
    extend (Z :. (1 :: Int) :. All :. All) .
    sumS . fromUnboxed (Z :. numPoint :. numPoint :. (1 :: Int)) $
    recon
  plotImageRepa (folderPath </> "Greens3DReconNorm_2.png") .
    ImageRepa 8 .
    computeS .
    reduceContrast 10 .
    extend (Z :. (1 :: Int) :. All :. All) .
    sumS . fromUnboxed (Z :. numPoint :. numPoint :. (1 :: Int)) $
    recon
