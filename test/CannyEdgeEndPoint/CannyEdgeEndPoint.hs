{-# LANGUAGE DataKinds #-}
module CannyEdgeEndPoint where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.ByteString           as BS (readFile)
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           GHC.Word
import           Image.IO
import           Image.IO
import           Image.Transform
import           Utils.Parallel
import           OpenCV                    as CV hiding (Z)
import           Data.Array.Unboxed        as AU
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array

main = do
  args <- getArgs
  let (numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFilePath:numIterationStr:writeSourceFlagStr:cutoffRadiusEndPointStr:cutoffRadiusStr:reversalFactorStr:inputImgPath:threshold1Str:threshold2Str:pixelDistStr:patchNormFlagStr:patchNormSizeStr:approximatedEigenValueStr:numThreadStr:_) =
        args
      numPoint = read numPointStr :: Int
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
      threshold1 = read threshold1Str :: Double
      threshold2 = read threshold2Str :: Double
      pixelDist = read pixelDistStr :: Int
      patchNormFlag = read patchNormFlagStr :: Bool
      patchNormSize = read patchNormSizeStr :: Int
      approximatedEigenValue = read approximatedEigenValueStr :: Double
      numThread = read numThreadStr :: Int
      folderPath = "output/test/CannyEdgeEndPoint"
  createDirectoryIfMissing True folderPath
  img <-
    (exceptError . coerceMat . imdecode ImreadUnchanged) <$>
    BS.readFile inputImgPath :: IO (Mat ('S '[ 'D, 'D]) 'D ('S GHC.Word.Word8))
  let edge =
        exceptError . canny threshold1 threshold2 (Just 7) CannyNormL2 $ img
      edgeRepa = R.map fromIntegral . toRepa $ edge
      (Z :. _ :. rows :. cols) = extent edgeRepa
      nonzeros = L.filter (> 0) . R.toList $ edgeRepa
      avg = L.sum nonzeros / (fromIntegral . L.length $ nonzeros)
      edgeRepaSparse =
        computeS .
        pad [numPoint, numPoint, 1] 0 .
        R.extend (Z :. (1 :: Int) :. All :. All) .
        increasePixelDistance pixelDist .
        R.slice
          (R.map
             (\x ->
                if x >= avg
                  then 255
                  else 0)
             edgeRepa) $
        (Z :. (0 :: Int) :. All :. All)
  plotImageRepa (folderPath </> "edge.png") . ImageRepa 8 . computeS $ edgeRepa
  plotImageRepa (folderPath </> "edge_sparse.png") . ImageRepa 8 $
    edgeRepaSparse
  edgeRepaSparse <-
    (\(ImageRepa _ img) -> img) <$>
    readImageRepa (folderPath </> "edge_sparse.png") False
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
  -- planReversal <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  -- (plan, pathNormMethod) <-
  --   makePatchNormFilter
  --     planReversal
  --     numPoint
  --     numPoint
  --     patchNormFlag
  --     patchNormSize
  plan <- makeR2Z2T0S0Plan emptyPlan arrR2Z2T0S0
  pathNormMethod <-
    if patchNormFlag
      then do
        let points =
              createIndex2D .
              L.map fst . L.filter (\(_, v) -> v /= 0) . AU.assocs $
              (AU.listArray ((0, 0), (numPoint - 1, numPoint - 1)) . R.toList $
               edgeRepaSparse :: AU.Array (Int, Int) Double)
            ys =
              pointCluster
                (connectionMatrixP
                   (ParallelParams numThread 1)
                   (2 * pixelDist)
                   points) $
              points
        M.zipWithM_
          (\i ->
             plotImageRepa (folderPath </> (printf "Cluster%03d.png" i)) .
             ImageRepa 8)
          [1 :: Int ..] .
          cluster2Array numPoint numPoint $
          ys
        return . PowerMethodConnection $ ys
      else return PowerMethodGlobal
  let endPointFilePath =
        folderPath </>
        (printf
           "EndPoint_%d_%d_%d_%d_%d_%d_%.2f_%.2f_%f.dat"
           numPoint
           (round thetaFreq :: Int)
           (round scaleFreq :: Int)
           (round maxScale :: Int)
           (round tao :: Int)
           cutoffRadiusEndPoint
           thetaSigma
           scaleSigma
           reversalFactor)
      edgeRepaSparseNormalized = normalizeValueRange (0, 1) edgeRepaSparse
      bias =
        R.traverse
          edgeRepaSparseNormalized
          (const
             (Z :. (L.length theta0Freqs) :. (L.length scale0Freqs) :. numPoint :.
              numPoint)) $ \f (Z :. _ :. _ :. i :. j) ->
          f (Z :. (0 :: Int) :. i :. j) :+ 0
      eigenVec =
        R.traverse
          edgeRepaSparseNormalized
          (const
             (Z :. (L.length thetaFreqs) :. (L.length scaleFreqs) :.
              (L.length theta0Freqs) :.
              (L.length scale0Freqs) :.
              numPoint :.
              numPoint)) $ \f (Z :. k :. l :. _ :. _ :. i :. j) ->
          if k == div (L.length thetaFreqs) 2 &&
             l == div (L.length scaleFreqs) 2 &&
             f (Z :. (0 :: Int) :. i :. j) > 0
            then 1 /
                 (fromIntegral $ (L.length theta0Freqs * L.length scale0Freqs)) :+
                 0
            else 0
  -- plotImageRepaComplex (folderPath </> "bias_0.png") .
  --   ImageRepa 8 .
  --   computeS .
  --   R.extend (Z :. (1 :: Int) :. All :. All) .
  --   R.sumS . R.sumS . rotate4D . rotate4D $
  --   bias
  endPointFlag <- doesFileExist endPointFilePath
  endPointR2Z2 <-
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
               completionFieldR2Z2 <-
                 computeS <$>
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
                   approximatedEigenValue
                   (computeS bias)
                   eigenVec
               -- writeRepaArray endPointFilePath completionFieldR2Z2
               return completionFieldR2Z2)
  let endPointR2Z2Biased = R.zipWith (*) endPointR2Z2 bias
      biasMag =
        R.sumS .
        R.sumS .
        rotate4D .
        rotate4D .
        R.map magnitude .
        r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
        endPointR2Z2Biased
  plotImageRepa (folderPath </> "EndPointBias.png") .
    ImageRepa 8 . computeS . R.extend (Z :. (1 :: Int) :. All :. All) $
    biasMag
  -- printf "%f %f\n" reversalFactor (R.sumAllS biasMag)
  -- let rotatedEndBias =
  --       computeUnboxedS $
  --       R.zipWith
  --         (+)
  --         (rotateBiasR2Z2T0S0 90 thetaFreqs endPointR2Z2Biased)
  --         (rotateBiasR2Z2T0S0 (-90) thetaFreqs endPointR2Z2Biased)
  -- powerMethodR2Z2T0S0Bias
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
  --   PowerMethodNone
  --   numIteration
  --   writeSourceFlag
  --   (printf
  --      "_%d_%d_%d_%d_%d_%d_%.2f_%.2f_%f"
  --      numPoint
  --      (round thetaFreq :: Int)
  --      (round scaleFreq :: Int)
  --      (round maxScale :: Int)
  --      (round tao :: Int)
  --      cutoffRadiusEndPoint
  --      thetaSigma
  --      scaleSigma
  --      reversalFactor)
  --   0.5
  --   rotatedEndBias
  --   eigenVec
