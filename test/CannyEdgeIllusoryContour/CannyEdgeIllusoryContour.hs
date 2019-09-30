{-# LANGUAGE DataKinds #-}
module CannyEdgeIllusoryContour where

import           Control.Monad             as M
import           Control.Monad.IO.Class
import           Data.Array.Repa           as R
import           Data.Array.Unboxed        as AU
import           Data.Binary
import           Data.ByteString           as BS
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.DomainChange
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           GHC.Word
import           Image.IO
import           Image.Transform
import           Linear.V2
import           OpenCV                    as CV hiding (Z)
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array
import           Utils.Parallel            hiding ((.|))
import           Utils.Time


main = do
  args <- getArgs
  let (numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:theta0FreqsStr:thetaFreqsStr:scale0FreqsStr:scaleFreqsStr:histFileName:histFileNameEndPoint:numIterationStr:numIterationEndPointStr:writeSourceFlagStr:cutoffRadiusEndPointStr:cutoffRadiusStr:reversalFactorStr:inputImgPath:threshold1Str:threshold2Str:writeSegmentsFlagStr:writeEndPointFlagStr:minSegLenStr:useFFTWWisdomFlagStr:fftwWisdomFileName:numThreadStr:_) =
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
      numIterationEndPoint = read numIterationEndPointStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      cutoffRadiusEndPoint = read cutoffRadiusEndPointStr :: Int
      cutoffRadius = read cutoffRadiusStr :: Int
      reversalFactor = read reversalFactorStr :: Double
      threshold1 = read threshold1Str :: Double
      threshold2 = read threshold2Str :: Double
      minSegLen = read minSegLenStr :: Int
      writeSegmentsFlag = read writeSegmentsFlagStr :: Bool
      writeEndPointFlag = read writeEndPointFlagStr :: Bool
      minimumPixelDist = 1 :: Int
      useFFTWWisdomFlag = read useFFTWWisdomFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/CannyEdgeIllusoryContour"
      histFilePath = folderPath </> histFileName
      histFilePathEndPoint = folderPath </> histFileNameEndPoint
      edgeFilePath = folderPath </> (takeBaseName inputImgPath L.++ "_edge.png")
      segmentsFilePath =
        folderPath </> (takeBaseName inputImgPath L.++ "_segments.dat")
      endPointFilePath =
        folderPath </>
        (printf
           "%s_EndPoint_%d_%d_%d_%d_%.2f_%f.dat"
           (takeBaseName inputImgPath)
           (numPoint * minimumPixelDist)
           (round thetaFreq :: Int)
           (round tao :: Int)
           cutoffRadiusEndPoint
           thetaSigma
           reversalFactor)
      fftwWisdomFilePath = folderPath </> fftwWisdomFileName
  createDirectoryIfMissing True folderPath
  copyFile inputImgPath (folderPath </> takeFileName inputImgPath)
  -- Find edges using Canny dege detector (Opencv)
  img <-
    (exceptError . coerceMat . imdecode ImreadUnchanged) <$>
    BS.readFile inputImgPath :: IO (Mat ('S '[ 'D, 'D]) 'D ('S GHC.Word.Word8))
  let [cols, rows] = miShape . matInfo $ img
      n = numPoint - cutoffRadiusEndPoint - 1 -- make sure that there is a zero wrap
      (w, h) =
        if rows > cols
          then ( n
               , round $ fromIntegral cols * fromIntegral n / fromIntegral rows)
          else ( round $ fromIntegral rows * fromIntegral n / fromIntegral cols
               , n)
      resizedImg =
        exceptError .
        resize
          (ResizeAbs . toSize $ V2 (fromIntegral w) (fromIntegral h))
          InterCubic $
        img
      edge =
        exceptError . canny threshold1 threshold2 (Just 3) CannyNormL2 $
        resizedImg
      edgeRepa =
        pad [numPoint, numPoint, 1] 0 .
        normalizeValueRange (0, 1) . R.map fromIntegral . toRepa $
        edge
  plotImageRepa edgeFilePath . ImageRepa 8 . computeS $ edgeRepa
  edgeRepa <- (\(ImageRepa _ img) -> img) <$> readImageRepa edgeFilePath False
  doesEndPointFileExist <- doesFileExist endPointFilePath
  endPointArray <-
    if doesEndPointFileExist && (not writeEndPointFlag)
      then do
        printCurrentTime "Read endpoint from file."
        readRepaArray endPointFilePath
      else do
        printCurrentTime "Start computing endpoint..."
        -- find segments for normalization:
        -- two adjacent non-zero-valued points are connected
        patchNormMethod <-
          if writeSegmentsFlag
            then do
              printCurrentTime "Start computing segments..."
              let nonzerorPoints =
                    createIndex2D .
                    L.map fst . L.filter (\(_, v) -> v /= 0) . AU.assocs $
                    (AU.listArray ((0, 0), (numPoint - 1, numPoint - 1)) .
                     R.toList $
                     edgeRepa :: AU.Array (Int, Int) Double)
                  segments =
                    L.map
                      (L.map
                         (\(a, b) ->
                            (a * minimumPixelDist, b * minimumPixelDist))) .
                    L.filter (\xs -> L.length xs >= minSegLen) .
                    pointCluster
                      (connectionMatrixP
                         (ParallelParams numThread 1)
                         1
                         nonzerorPoints) $
                    nonzerorPoints
              encodeFile segmentsFilePath segments
              removePathForcibly (folderPath </> "segments")
              createDirectoryIfMissing True (folderPath </> "segments")
              M.zipWithM_
                (\i ->
                   plotImageRepa
                     (folderPath </> "segments" </> (printf "Cluster%03d.png" i)) .
                   ImageRepa 8)
                [1 :: Int ..] .
                cluster2Array
                  (numPoint * minimumPixelDist)
                  (numPoint * minimumPixelDist) $
                segments
              printCurrentTime "Done computing segments."
              return . PowerMethodConnection $ segments
            else do
              printCurrentTime "Read segments from files."
              PowerMethodConnection <$> decodeFile segmentsFilePath
        -- Compute the Green's function
        let numPointEndPoint = minimumPixelDist * numPoint
            maxScaleEndPoint = 1.00000000001
        flag <- doesFileExist histFilePathEndPoint
        radialArr <-
          if flag
            then do
              printCurrentTime "Read endpoint filter histogram from files."
              R.map magnitude . getNormalizedHistogramArr <$>
                decodeFile histFilePathEndPoint
            else do
              printCurrentTime
                "Couldn't find a Green's function data. Start simulation..."
              solveMonteCarloR2Z2T0S0Radial
                numThread
                numTrail
                maxTrail
                numPointEndPoint
                numPointEndPoint
                thetaSigma
                0.0
                maxScaleEndPoint
                tao
                theta0Freqs
                thetaFreqs
                [0]
                [0]
                histFilePathEndPoint
                (emptyHistogram
                   [ (round . sqrt . fromIntegral $
                      2 * (div numPointEndPoint 2) ^ 2)
                   , 1
                   , L.length theta0Freqs
                   , 1
                   , L.length thetaFreqs
                   ]
                   0)
        arrR2Z2T0S0 <-
          computeUnboxedP $
          computeR2Z2T0S0ArrayRadial
            (pinwheelHollowNonzeronCenter 12)
            (cutoff cutoffRadiusEndPoint radialArr)
            numPointEndPoint
            numPointEndPoint
            1
            maxScaleEndPoint
            thetaFreqs
            [0]
            theta0Freqs
            [0]
        plan <-
          makeR2Z2T0S0Plan
            emptyPlan
            useFFTWWisdomFlag
            fftwWisdomFilePath
            arrR2Z2T0S0
        -- Compute initial eigenvector and bias
        -- increase the distance between pixels to aovid aliasing
        let resizedEdgeRepa =
              R.traverse
                edgeRepa
                (const (Z :. numPointEndPoint :. numPointEndPoint)) $ \f (Z :. i :. j) ->
                if mod i minimumPixelDist == 0 && mod j minimumPixelDist == 0
                  then f (Z :. (0 :: Int) :. (div i minimumPixelDist) :.
                          (div j minimumPixelDist))
                  else 0
            bias =
              computeBiasR2T0S0FromRepa
                numPointEndPoint
                numPointEndPoint
                (L.length theta0Freqs)
                1
                resizedEdgeRepa
            eigenVec =
              computeInitialEigenVectorR2T0S0FromRepa
                numPointEndPoint
                numPointEndPoint
                (L.length theta0Freqs)
                1
                (L.length thetaFreqs)
                1
                resizedEdgeRepa
        plotImageRepa
          (folderPath </> takeBaseName inputImgPath L.++ "_resizedEdge.png") .
          ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
          resizedEdgeRepa
        endPointSource <-
          computeS . R.zipWith (*) bias <$>
          powerMethodR2Z2T0S0Reversal
            plan
            folderPath
            numPointEndPoint
            numPointEndPoint
            numOrientation
            thetaFreqs
            theta0Freqs
            1
            [0]
            [0]
            0
            arrR2Z2T0S0
            patchNormMethod
            numIterationEndPoint
            writeSourceFlag
            (printf
               "_%d_%d_%d_%d_%.2f_%f_%s_EndPoint"
               (numPoint * minimumPixelDist)
               (round thetaFreq :: Int)
               (round tao :: Int)
               cutoffRadiusEndPoint
               thetaSigma
               reversalFactor
               (takeBaseName inputImgPath))
            0.5
            reversalFactor
            bias
            eigenVec
        writeRepaArray endPointFilePath endPointSource
        return endPointSource
  print "done"

