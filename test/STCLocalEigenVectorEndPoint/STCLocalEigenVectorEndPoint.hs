module STCLocalEigenVectorEndPoint where

import           Control.Monad             as M
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.MonteCarlo
import           FokkerPlanck.Pinwheel
import           Image.IO
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Types
import           Utils.Array
import           Utils.Parallel
import           Filter.Pinwheel

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:thetaFreqsStr:scaleFreqsStr:histFileName:numIterationStr:writeSourceFlagStr:cutoffRadiusStr:patchNormFlagStr:patchNormSizeStr:shape2DStr:useFFTWWisdomFlagStr:fftwWisdomFileName:numThreadStr:_) <-
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
      thetaFreq = read thetaFreqsStr :: Double
      thetaFreqs = [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Double
      scaleFreqs = [-scaleFreq .. scaleFreq]
      numIteration = read numIterationStr :: Int
      writeSourceFlag = read writeSourceFlagStr :: Bool
      cutoffRadius = read cutoffRadiusStr :: Int
      patchNormFlag = read patchNormFlagStr :: Bool
      patchNormSize = read patchNormSizeStr :: Int
      shape2D@(Points _ minDist _) = read shape2DStr :: Points Shape2D
      useFFTWWisdomFlag = read useFFTWWisdomFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCLocalEigenVectorEndPoint"
      histFilePath = folderPath </> histFileName
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
          thetaFreqs
          thetaFreqs
          scaleFreqs
          scaleFreqs
          histFilePath
          (emptyHistogram
             [ (round . sqrt . fromIntegral $ 2 * (div numPoint 2) ^ 2)
             , L.length scaleFreqs
             , L.length thetaFreqs
             , L.length scaleFreqs
             , L.length thetaFreqs
             ]
             0)
  let localEigenVector =
        computeLocalEigenVector
          (ParallelParams numThread 0)
          (pinwheelHollow 2)
          (cutoff cutoffRadius radialArr)
          numPoint
          numPoint
          maxScale
          thetaFreqs
          scaleFreqs
      pinwheelArr =
        -- computeS .
        R.traverse2
          (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
          (fromListUnboxed (Z :. (L.length scaleFreqs)) scaleFreqs)
          (\(Z :. numThetaFreq) (Z :. numScaleFreq) ->
             (Z :. numThetaFreq :. numScaleFreq :. numPoint :. numPoint)) $ \ft fs (Z :. k :. l :. i :. j) ->
          let r =
                sqrt . fromIntegral $
                (i - center numPoint) ^ 2 + (j - center numPoint) ^ 2
           in if r < fromIntegral cutoffRadius
                then 0
                else pinwheelHollow
                       1
                       (ft (Z :. k))
                       (fs (Z :. l))
                       maxScale
                       0
                       (i - center numPoint)
                       (j - center numPoint)
      -- localEigenVector = (computeS $
      --                     R.traverse2
      --                       localEigenVector'
      --                       (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
      --                       const $ \f fFreq idx@(Z :. k :. _ :. _ :. _) ->
      --                       f idx * exp (0 :+ fFreq (Z :. k) * pi))
      pointSet = makeShape2D shape2D
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
      bias = computeBiasR2T0S0 numPoint numPoint thetaFreqs scaleFreqs xs
  arrR2Z2T0S0 <-
    computeUnboxedP $
    computeR2Z2T0S0ArrayRadial
      (pinwheelHollow 4)
      (cutoff cutoffRadius radialArr)
      numPoint
      numPoint
      1
      maxScale
      thetaFreqs
      scaleFreqs
      thetaFreqs
      scaleFreqs
  plan' <-
    make4DPlan emptyPlan useFFTWWisdomFlag fftwWisdomFilePath localEigenVector
  plan <-
    makeR2Z2T0S0Plan plan' useFFTWWisdomFlag fftwWisdomFilePath arrR2Z2T0S0
  plotImageRepa (folderPath </> "Shape.png") . ImageRepa 8 $ shapeArr
  filterF <- dft4D plan . computeS . makeFilter2D $ localEigenVector
  initialVec <-
    R.zipWith (*) bias <$>
    convolve4D
      plan
      filterF
      (computeInitialDistributionR2T0S0'
         numPoint
         numPoint
         thetaFreqs
         scaleFreqs
         maxScale
         xs)
  -- eigenVec <-
  --   powerMethodR2Z2
  --     plan
  --     folderPath
  --     numPoint
  --     numPoint
  --     numOrientation
  --     thetaFreqs
  --     numScale
  --     scaleFreqs
  --     maxScale
  --     pinwheelArr
  --     -- localEigenVector
  --     PowerMethodGlobal
  --     numIteration
  --     writeSourceFlag
  --     (printf
  --        "_%d_%d_%d_%d_%d_%d_%.2f_%.2f"
  --        numPoint
  --        (round thetaFreq :: Int)
  --        (round scaleFreq :: Int)
  --        (round maxScale :: Int)
  --        (round tao :: Int)
  --        cutoffRadius
  --        thetaSigma
  --        scaleSigma)
  --     bias
  --     initialVec
  --     -- bias
  --     -- (R.traverse bias id $ \f idx@(Z :. k :. l :. _ :. _) ->
  --     --    if k == div (L.length thetaFreqs) 2 && l == div (L.length scaleFreqs) 2
  --     --      then f idx
  --     --      else 0)
  pathNormMethod <-
    if patchNormFlag
      then do
        let points =
              createIndex2D .
              L.map (\(i, j) -> (i + div numPoint 2, j + div numPoint 2)) .
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
             plotImageRepa (folderPath </> (printf "Cluster%03d.png" i)) .
             ImageRepa 8)
          [1 :: Int ..] .
          cluster2Array numPoint numPoint $
          ys
        return . PowerMethodConnection $ ys
      else return PowerMethodGlobal
  sourceR2Z2 <-
    powerMethodR2Z2T0S0Reversal
      plan
      folderPath
      numPoint
      numPoint
      numOrientation
      thetaFreqs
      thetaFreqs
      numScale
      scaleFreqs
      scaleFreqs
      maxScale
      arrR2Z2T0S0
      pathNormMethod
      0
      writeSourceFlag
      (printf
         "_%d_%d_%d_%d_%d_%d_%.2f_%.2f_Final"
         numPoint
         (round thetaFreq :: Int)
         (round scaleFreq :: Int)
         (round maxScale :: Int)
         (round tao :: Int)
         cutoffRadius
         thetaSigma
         scaleSigma)
      0.0
      0.0
      bias
      -- (computeInitialEigenVectorR2T0S0
      --    numPoint
      --    numPoint
      --    thetaFreqs
      --    scaleFreqs
      --    thetaFreqs
      --    scaleFreqs
      --    xs)
      (R.traverse
         -- eigenVec
         initialVec
         (\(Z :. numThetaFreq :. numScaleFreq :. x :. y) ->
            (Z :. numThetaFreq :. numScaleFreq :. numThetaFreq :. numScaleFreq :.
             x :.
             y)) $ \f idx@(Z :. tf :. sf :. t0f :. s0f :. i :. j) ->
         f (Z :. tf :. sf :. i :. j))
  powerMethodR2Z2
    plan
    folderPath
    numPoint
    numPoint
    numOrientation
    thetaFreqs
    numScale
    scaleFreqs
    maxScale
    localEigenVector
    PowerMethodGlobal
    numIteration
    writeSourceFlag
    (printf
       "_%d_%d_%d_%d_%d_%d_%.2f_%.2f_TEST"
       numPoint
       (round thetaFreq :: Int)
       (round scaleFreq :: Int)
       (round maxScale :: Int)
       (round tao :: Int)
       cutoffRadius
       thetaSigma
       scaleSigma)
    bias
    (R.zipWith (*) bias sourceR2Z2)
