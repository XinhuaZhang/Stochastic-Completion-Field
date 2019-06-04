module STCR2Z2T0S0EndPoint where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           Image.Transform           (normalizeValueRange)
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           System.Random
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:cutoffRadiusEndPointStr:cutoffRadiusStr:reversalFactorStr:cStr:numThreadStr:_) <-
    getArgs
  print args
  let numPoint = read numPointStr :: Int
      numOrientation = read numOrientationStr :: Int
      numScale = read numScaleStr :: Int
      thetaSigma = read thetaSigmaStr :: Double
      scaleSigma = read scaleSigmaStr :: Double
      maxScale = read maxScaleStr :: Double
      tao = read taoStr :: Double
      numTrail = read numTrailStr :: Int
      maxTrail = read maxTrailStr :: Int
      theta0Freq = read theta0FreqsStr :: Double
      theta0Freqs = [-theta0Freq .. theta0Freq]
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scale0Freq = read scale0FreqsStr :: Double
      scaleFreq = read scaleFreqsStr :: Double
      scale0Freqs = [-scale0Freq .. scale0Freq]
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      cutoffRadiusEndPoint = read cutoffRadiusEndPointStr :: Int
      cutoffRadius = read cutoffRadiusStr :: Int
      reversalFactor = read reversalFactorStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCR2Z2T0S0EndPoint"
      a = 20 :: Int
      b = 5 :: Int
      c = read cStr :: Int
      endPointFilePath =
        folderPath </>
        (printf
           "EndPoint_%d_%d_%d_%d_%d_%d_%.2f_%.2f_%d_%d_%d_%f.dat"
           numPoint
           (round thetaFreq :: Int)
           (round scaleFreq :: Int)
           (round maxScale :: Int)
           (round tao :: Int)
           cutoffRadiusEndPoint
           thetaSigma
           scaleSigma
           a
           b
           c
           reversalFactor)
  createDirectoryIfMissing True folderPath
  flag <- doesFileExist histFilePath
  radialArr <-
    if flag
      then R.map magnitude . getNormalizedHistogramArr <$>
           decodeFile histFilePath
      else do
        putStrLn "Couldn't find a Green's function data. Start simulation..."
        solveMonteCarloR2Z2T0S0Radial
          numThread
          numTrail
          maxTrail
          numPoint
          numPoint
          thetaSigma
          scaleSigma
          maxScale
          tao
          theta0Freqs
          thetaFreqs
          scale0Freqs
          scaleFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scale0Freqs
             , L.length theta0Freqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (cutoff cutoffRadius radialArr)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
      -- a' = round $ (fromIntegral a) * (sqrt 2) / 2
      -- b' = round $ (fromIntegral b) * (sqrt 2) / 2
      -- c' = round $ (fromIntegral c) * (sqrt 2) / 2
      -- xs =
      --   ([R2S1RPPoint (i, i, 0, 0) | i <- [a',a' + b' .. c']] L.++
      --    [R2S1RPPoint (i, -i, 0, 0) | i <- [-a',-(a' + b') .. -c']] L.++
      --    [R2S1RPPoint (i, i, 0, 0) | i <- [-a',-(a' + b') .. -c']] L.++
      --    [R2S1RPPoint (i, -i, 0, 0) | i <- [a',a' + b' .. c']] L.++
      --    [R2S1RPPoint (i, 0, 0, 0) | i <- [a,a + b .. c]] L.++
      --    [R2S1RPPoint (0, i, 0, 0) | i <- [-a,-(a + b) .. -c]] L.++
      --    [R2S1RPPoint (i, 0, 0, 0) | i <- [-a,-(a + b) .. -c]] L.++
      --    [R2S1RPPoint (0, i, 0, 0) | i <- [a,a + b .. c]])
      -- xs = [R2S1RPPoint (i, 0, 0, 0) | i <- [a,a + b .. c]] L.++ [R2S1RPPoint (0, i, 0, 0) | i <- [b,2 * b .. c]] L.++ [R2S1RPPoint (0, -i, 0, 0) | i <- [b,2 * b .. c]]
      -- r = 50 :: Double
      -- numTheta = 26 :: Double
      -- theta = 2 * pi - pi / 3
      -- deltaTheta = theta / numTheta
      -- xs =
      --   (L.map
      --      (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
      --      [ (r * cos (k * deltaTheta) + 0, r * sin (k * deltaTheta) + 0)
      --      | k <- [0 .. numTheta]
      --      ]) L.++
      --   [R2S1RPPoint (i, 0, 0, 0) | i <- [0,10 .. round r - 10]] L.++
      --   (L.map
      --      (\(i, j) -> R2S1RPPoint (round i, round j, 0, 1))
      --      [(r' * cos theta + 0, r' * sin theta + 0) | r' <- [10,20 .. r - 10]])
  xs <-
    M.foldM
      (\((R2S1RPPoint (z0, 0, 0, 0)):ys) _ -> do
         z <- randomRIO (4, 10) :: IO Double
         if z0 + round z <= 25
           then return $
                (R2S1RPPoint (z0 + round z, 0, 0, 0)) :
                ((R2S1RPPoint (z0, 0, 0, 0)) : ys)
           else return ((R2S1RPPoint (z0, 0, 0, 0)) : ys))
      [R2S1RPPoint (-25, 0, 0, 0)]
      [1 .. 7]
  let bias = computeBiasR2T0S0 numPoint numPoint theta0Freqs scale0Freqs xs
      eigenVec =
        computeInitialEigenVectorR2T0S0
          numPoint
          numPoint
          theta0Freqs
          scale0Freqs
          thetaFreqs
          scaleFreqs
          xs
  endPointFlag <- doesFileExist endPointFilePath
  completionFieldR2Z2'' <-
    if endPointFlag
      then readRepaArray endPointFilePath
      else (do putStrLn "Couldn't find the endpoint data. Start computing..."
               arrR2Z2T0S0EndPoint <-
                 computeUnboxedP $
                 computeR2Z2T0S0ArrayRadial
                   (cutoff cutoffRadiusEndPoint radialArr)
                   numPoint
                   numPoint
                   1
                   maxScale
                   thetaFreqs
                   scaleFreqs
                   theta0Freqs
                   scale0Freqs
               (plan', pathNormFilterF) <-
                 makePatchNormFilter plan numPoint numPoint 15
               completionFieldR2Z2' <-
                 computeS <$>
                 powerMethodR2Z2T0S0Reversal
                   plan'
                   folderPath
                   numPoint
                   numPoint
                   numOrientation
                   thetaFreqs
                   theta0Freqs
                   numScale
                   scaleFreqs
                   scale0Freqs
                   arrR2Z2T0S0EndPoint
                   pathNormFilterF
                   numIteration
                   writeSourceFlag
                   (printf
                      "_%d_%d_%d_%d_%d_%d_%.2f_%.2f_%f_EndPoint"
                      numPoint
                      (round thetaFreq :: Int)
                      (round scaleFreq :: Int)
                      (round maxScale :: Int)
                      (round tao :: Int)
                      cutoffRadiusEndPoint
                      thetaSigma
                      scaleSigma
                      reversalFactor)
                   0.5
                   reversalFactor
                   bias
                   eigenVec
               -- writeRepaArray endPointFilePath completionFieldR2Z2'
               return completionFieldR2Z2')
  let completionFieldR2Z2 = R.zipWith (*) completionFieldR2Z2'' bias
      endpointBias = completionFieldR2Z2
        -- rotateBiasR2Z2T0S0 0 theta0Freqs . R.traverse completionFieldR2Z2 id $ \f idx@(Z :. _ :. _ :. i :. j) ->
        --   if (sqrt . fromIntegral $
        --       (i - div numPoint 2) ^ 2 + (j - div numPoint 2) ^ 2) >
        --      35
        --     then 0
        --     else f idx
        -- R.traverse
          -- (R.zipWith
          --    (+)
          --    (rotateBiasR2Z2T0S0 90 theta0Freqs completionFieldR2Z2)
          --    (rotateBiasR2Z2T0S0 (-90) theta0Freqs completionFieldR2Z2))
        --   id $ \f idx@(Z :. _ :. _ :. i :. j) ->
        --   if (sqrt . fromIntegral $
        --       (i - div numPoint 2) ^ 2 + (j - div numPoint 2) ^ 2) >
        --      35
        --     then 0
        --     else f idx
      biasMag =
        R.sumS .
        R.sumS .
        rotate4D .
        rotate4D .
        R.map magnitude .
        r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
        endpointBias
  plotImageRepa (folderPath </> "EndPointBias.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    biasMag
  printf "%f %f\n" reversalFactor (R.sumAllS biasMag)
  -- powerMethodR2Z2T0S0BiasReversal
  --   plan
  --   folderPath
  --   numPoint
  --   numPoint
  --   numOrientation
  --   thetaFreqs
  --   theta0Freqs
  --   numScale
  --   scaleFreqs
  --   scale0Freqs
  --   arrR2Z2T0S0
  --   -- numIteration
  --   10
  --   writeSourceFlag
  --   -- True
  --   (printf
  --      "_%d_%d_%d_%d_%d_%d_%.2f_%.2f"
  --      numPoint
  --      (round thetaFreq :: Int)
  --      (round scaleFreq :: Int)
  --      (round maxScale :: Int)
  --      (round tao :: Int)
  --      cutoffRadius
  --      thetaSigma
  --      scaleSigma)
  --   1
  --   (computeS endpointBias)
  --   (R.fromFunction
  --      (Z :. (L.length thetaFreqs) :. (L.length scaleFreqs) :.
  --       (L.length theta0Freqs) :.
  --       (L.length scale0Freqs) :.
  --       numPoint :.
  --       numPoint) $ \(Z :. k :. l :. _ :. _ :. i :. j) ->
  --      if k == div (L.length thetaFreqs) 2 && l == div (L.length scaleFreqs) 2
  --        then 1 / (fromIntegral $ (L.length theta0Freqs * L.length scale0Freqs)) :+
  --             0
  --        else 0)
