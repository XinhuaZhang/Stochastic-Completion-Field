module STCLocalEigenVector where

import           Control.Monad             as M
import           Control.Monad.Parallel             as MP
import           Data.Array.Repa           as R
import           Data.Binary               (decodeFile)
import           Data.Complex
import           Data.List                 as L
import           DFT.Plan
import           FokkerPlanck.DomainChange
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

{-# INLINE timeReversal #-}
timeReversal :: (R.Source r e) => Array r DIM4 e -> Array D DIM4 e
timeReversal arr =
  let (Z :. nf :. _ :. _ :. _) = extent arr
      n = div nf 2
   in R.backpermute
        (extent arr)
        (\(Z :. k :. l :. i :. j) ->
           let x = nf - n
            in if k >= x
                 then (Z :. k - x :. l  :. i :. j)
                 else (Z :. k + n :. l :. i :. j))
        arr

main = do
  args@(numPointStr:numOrientationStr:numScaleStr:thetaSigmaStr:scaleSigmaStr:maxScaleStr:taoStr:numTrailStr:maxTrailStr:thetaFreqsStr:scaleFreqsStr:histFileName:cutoffRadiusStr:initDistStr:useFFTWWisdomFlagStr:fftwWisdomFileName:numThreadStr:_) <-
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
      cutoffRadius = read cutoffRadiusStr :: Int
      initDist = read initDistStr :: [R2S1RPPoint]
      useFFTWWisdomFlag = read useFFTWWisdomFlagStr :: Bool
      numThread = read numThreadStr :: Int
      folderPath = "output/test/STCLocalEigenVector"
      histFilePath = folderPath </> histFileName
      fftwWisdomFilePath = folderPath </> fftwWisdomFileName
      sourceDist = L.take 1 initDist
      sinkDist = L.drop 1 initDist
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
      localEigenVectorSink =
        computeLocalEigenVectorSink
          (ParallelParams numThread 0)
          (pinwheelHollow 2)
          (cutoff cutoffRadius radialArr)
          numPoint
          numPoint
          maxScale
          thetaFreqs
          scaleFreqs
  plan <-
    make4DPlan emptyPlan useFFTWWisdomFlag fftwWisdomFilePath localEigenVector
  sourceDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      thetaFreqs
      scaleFreqs
      maxScale
      sourceDist
  sinkDistArr <-
    computeInitialDistributionR2T0S0
      plan
      numPoint
      numPoint
      thetaFreqs
      scaleFreqs
      maxScale
      sinkDist
  filterF <- dft4D plan . computeS . makeFilter4D $ localEigenVector
  filterSinkF <- dft4D plan . computeS . makeFilter4D $ localEigenVector
  -- Source
  sourceArr <- convolve4D plan filterF sourceDistArr
  sourceR2 <-
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    sourceArr
  plotImageRepaComplex (folderPath </> "Source.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    sourceR2
  -- Sink
  sinkArr <- convolve4D plan filterSinkF sinkDistArr
  -- let sinkArr =
  --       R.traverse2
  --         sinkArr'
  --         (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
  --         const $ \f fFreq idx@(Z :. k :. l :. i :. j) ->
  --         f idx * ((exp (0 :+ (-1) * (fFreq (Z :. k)) * pi)))
  sinkR2 <-
    R.sumP .
    R.sumS .
    rotate4D .
    rotate4D . r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs $
    sinkArr
  plotImageRepaComplex (folderPath </> "Sink.png") .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    sinkR2
  -- Completion
  -- plotImageRepa (folderPath </> "Completion.png") .
  --   ImageRepa 8 .
  --   computeS .
  --   extend (Z :. (1 :: Int) :. All :. All) .
  --   R.sumS . R.sumS . rotate4D . rotate4D . R.map magnitude $
  --   (timeReversal $  r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs sinkArr) *^
  --   (r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs sourceArr)
  -- completionFieldR2 plan folderPath "" sourceR2 sinkR2
  completionFieldR2Z2
    plan
    folderPath
    ""
    numOrientation
    thetaFreqs
    numScale
    scaleFreqs
    sourceArr
    sinkArr
    -- (computeS $
    --  R.traverse2
    --    sinkArr
    --    (fromListUnboxed (Z :. (L.length thetaFreqs)) thetaFreqs)
    --    const $ \f fFreq idx@(Z :. k :. l :. i :. j) ->
    --    f idx * ((exp (0 :+ (-1) * (fFreq (Z :. k)) * pi))))
    -- (computeS . R.map (\x -> let (m,p) = polar x
    --                          in mkPolar m (p + pi)) $ sinkArr)
    -- (computeS $ flip4D sinkArr)
    -- (computeS $ timeReverse4D thetaFreqs sinkArr)
  -- MP.mapM_
  --   (\(tf, sf) ->
  --      plotImageRepaComplex (folderPath </> printf "Source_%d_%d.png" tf sf) .
  --      ImageRepa 8 .
  --      computeS .
  --      extend (Z :. (1 :: Int) :. All :. All) .
  --      R.slice
  --        (r2z2Tor2s1rp numOrientation thetaFreqs numScale scaleFreqs sourceArr) $
  --      (Z :. tf :. sf :. All :. All))
  --   [ (tf, sf)
  --   | tf <- [0 .. numOrientation - 1]
  --   , sf <- [0 .. numScale - 1]
  --   ]
    -- [ (tf, sf)
    -- | tf <- [0 .. L.length thetaFreqs - 1]
    -- , sf <- [0 .. L.length scaleFreqs - 1]
    -- ]
