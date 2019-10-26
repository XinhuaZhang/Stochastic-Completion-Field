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
import           Utils.Parallel

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:cutoffRadiusEndPointStr:cutoffRadiusStr:reversalFactorStr:cStr:patchNormFlagStr:patchNormSizeStr:approximatedEigenValueStr:shape2DStr:segmentsFilePath:segIdxStr:useFFTWWisdomFlagStr:fftwWisdomFileName:numThreadStr:_) <-
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
      patchNormFlag = read patchNormFlagStr :: Bool
      patchNormSize = read patchNormSizeStr :: Int
      approximatedEigenValue = read approximatedEigenValueStr :: Double
      shape2D@(Points _ minDist _) = read shape2DStr :: Points Shape2D
      segIdx = read segIdxStr :: Int
      useFFTWWisdomFlag = read useFFTWWisdomFlagStr :: Bool
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
      fftwWisdomFilePath = folderPath </> fftwWisdomFileName
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
      (PinwheelHollow0 10)
      (cutoff cutoffRadius radialArr)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      theta0Freqs
      scale0Freqs
  plan <-
    makeR2Z2T0S0Plan emptyPlan useFFTWWisdomFlag fftwWisdomFilePath arrR2Z2T0S0
  -- (plan, pathNormMethod) <-
  --   makePatchNormFilter plan' numPoint numPoint patchNormFlag patchNormSize
      -- minDist = 8
      -- kanizsaTriangle1 =
      --   makeShape2D $ Points (-30, -30) minDist (Corner 30 60 80) --(PacMan 30 60 50)  --(Ehrenstein 8 15 40)  -- (IncompleteCircle 0 60 50 ) -- (TJunction 45 50) -- (PacMan 0 60 100 ) -- (Corner 0 60 100) --(KanizsaTriangle1 0 480 160 80)
  -- ys <-
  --   (\aa -> aa L.!! segIdx) <$> decodeFile segmentsFilePath :: IO [(Int, Int)]
  let pointSet = makeShape2D shape2D
      shapeArr =
        getShape2DRepaArray
          numPoint
          numPoint
          (L.map
             (\(x, y) ->
                (x + fromIntegral numPoint / 2, y + fromIntegral numPoint / 2))
             pointSet)
      xs =
        L.map (\(x, y) -> R2S1RPPoint (x, y, 0, 1)) . getShape2DIndexList $
        pointSet
      -- (xAvg, yAvg) =
      --   (\(as, bs) ->
      --      ( round $
      --        (fromIntegral . L.sum $ as) / (fromIntegral . L.length $ as)
      --      , round $
      --        (fromIntegral . L.sum $ bs) / (fromIntegral . L.length $ bs))) .
      --   L.unzip $
      --   ys
      -- centeredYs = L.map (\(x, y) -> (x - xAvg, y - yAvg)) ys
      -- xs = L.map (\(x, y) -> R2S1RPPoint (x, y, 0, 1)) centeredYs
      -- shapeArr =
      --   L.head . cluster2Array numPoint numPoint $
      --   [L.map (\(x, y) -> (x + div numPoint 2, y + div numPoint 2)) centeredYs]
      -- pointSet = L.map (\(x, y) -> (fromIntegral x, fromIntegral y)) centeredYs
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
  plotImageRepa (folderPath </> "Shape.png") . ImageRepa 8 $ shapeArr
  endPointFlag <- doesFileExist endPointFilePath
  completionFieldR2Z2'' <-
    if endPointFlag
      then readRepaArray endPointFilePath
      else (do putStrLn "Couldn't find the endpoint data. Start computing..."
               arrR2Z2T0S0EndPoint <-
                 computeUnboxedP $
                 computeR2Z2T0S0ArrayRadial
                   -- pinwheel
                   -- (pinwheelHollowNonzeronCenter 16) 
                   (PinwheelHollow0 4)
                   (cutoff cutoffRadiusEndPoint radialArr)
                   numPoint
                   numPoint
                   1
                   maxScale
                   thetaFreqs
                   scaleFreqs
                   theta0Freqs
                   scale0Freqs
               pathNormMethod <-
                 if patchNormFlag
                   then do
                     let points =
                           createIndex2D .
                           L.map
                             (\(i, j) ->
                                (i + div numPoint 2, j + div numPoint 2)) .
                           getShape2DIndexList $
                           pointSet
                         ys =
                           pointCluster
                             (connectionMatrixP
                                (ParallelParams numThread 1)
                                (minDist + 1)
                                points) $
                           points
                     M.zipWithM_
                       (\i ->
                          plotImageRepa
                            (folderPath </> (printf "Cluster%03d.png" i)) .
                          ImageRepa 8)
                       [1 :: Int ..] .
                       cluster2Array numPoint numPoint $
                       ys
                     return . PowerMethodConnection $ ys
                   else return PowerMethodGlobal
               (R.foldAllP max 0 . R.map magnitude $ arrR2Z2T0S0EndPoint) >>=
                 print
               completionFieldR2Z2' <-
                 powerMethodR2Z2T0S0Reversal
                   plan
                   folderPath
                   numPoint
                   numPoint
                   numOrientation
                   thetaFreqs
                   theta0Freqs
                   numScale
                   scaleFreqs
                   scale0Freqs
                   maxScale
                   arrR2Z2T0S0EndPoint
                   pathNormMethod
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
      endpointBias = rotateBiasR2Z2T0S0 180 theta0Freqs completionFieldR2Z2
      -- endpointBias =
      --   R.zipWith
      --     (+)
      --     (rotateBiasR2Z2T0S0 90 theta0Freqs completionFieldR2Z2)
      --     (rotateBiasR2Z2T0S0 (-90) theta0Freqs completionFieldR2Z2)
                    -- rotateBiasR2Z2T0S0 0 theta0Freqs . R.traverse completionFieldR2Z2 id $ \f idx@(Z :. _ :. _ :. i :. j) ->
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
  --   0.1
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
