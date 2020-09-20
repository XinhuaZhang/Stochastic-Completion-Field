{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Strict       #-}
module STC.PowerMethod where

import           Control.Monad                  as M
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX as A
import           Data.Array.IArray              as IA
import           Data.Array.Repa                as R
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.FourierSeries
import           FourierMethod.BlockMatrixAcc
import           FourierMethod.FourierSeries2D
import           Graphics.Gnuplot.Simple
import           Image.IO
import           STC.CompletionField
import           STC.Convolution
import           STC.DFTArray
import           STC.Multiplication
import           STC.Plan
import           STC.Reversal
import           STC.Utils
import           System.FilePath                ((</>))
import           Text.Printf
import           Utils.Array
import           Utils.List
import           Utils.Parallel
import           Utils.Time
import FourierPinwheel as FP
import           Foreign.CUDA.Driver            as CUDA

powerMethod ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> DFTArray
  -> IO DFTArray
powerMethod _ _ _ _ _ _ 0 !arr = return arr
powerMethod !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !numStep !input@(DFTArray rows cols _ _ _) = do
  printCurrentTime (show numStep)
  convolvedArr <-
    STC.Convolution.convolve Source plan coefficients harmonicsArray input
  let !biasedConvolvedArr = parMapDFTArray (VS.zipWith (*) bias) convolvedArr
      !maxMag =
        L.sum .
        parMap rdeepseq (VS.sum . VS.map (\x -> (magnitude x)^2)) . getDFTArrayVector $
        biasedConvolvedArr
      !normalizedBiasedConvolvedArr =
        parMapDFTArray (VS.map (/ (maxMag :+ 0))) biasedConvolvedArr
  when writeFlag .
    plotDFTArrayPower
      (folderPath </> (printf "Source_%03d.png" numStep))
      cols
      rows $
    convolvedArr
  powerMethod
    plan
    folderPath
    writeFlag
    coefficients
    harmonicsArray
    bias
    (numStep - 1)
    normalizedBiasedConvolvedArr

{-# INLINE computeContour #-}
computeContour ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> String
  -> DFTArray
  -> IO DFTArray
computeContour !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !numIteration !suffix !input@(DFTArray rows cols _ _ _) = do
  eigenVector <-
    powerMethod
      plan
      folderPath
      writeFlag
      coefficients
      harmonicsArray
      bias
      numIteration
      input
  source <- STC.Convolution.convolve Source plan coefficients harmonicsArray eigenVector
  sink <- STC.Convolution.convolve Sink plan coefficients harmonicsArray eigenVector
  completion <- completionField plan source sink
  plotDFTArrayPower (folderPath </> "SinkPower.png") rows cols sink
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower_%s.png" suffix))
    rows
    cols
    completion
  return completion

powerMethod' ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int -> [[(Complex Double)]]
  -> DFTArray
  -> IO DFTArray
powerMethod' _ _ _ _ _ _ 0 _ !arr = return arr
powerMethod' !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !numStep thetaRHarmonics !input@(DFTArray rows cols thetaFreqs rFreqs _) = do
  printCurrentTime (show numStep)
  convolvedArr <- convolve' Source plan coefficients harmonicsArray input
  let -- sigma = 15
      -- gaussian =
      --   L.map (\x -> (exp (-(x ^ 2) / (2 * sigma ^ 2))) :+ 0) thetaFreqs
      convolvedArr' =
        fromUnboxed (Z :. 72 :. cols :. rows) .
        VS.convert .
        VS.concat .
        computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $
        convolvedArr
      -- convolvedArr1 =
      --   DFTArray rows cols thetaFreqs rFreqs .
      --   L.zipWith (\g vec -> VS.map (* g) vec) gaussian . getDFTArrayVector $
      --   convolvedArr
  plotThetaDimension folderPath (printf "FourierSeries_%d_" numStep) (29, 12) .
    R.map magnitude $
    convolvedArr'
  let !biasedConvolvedArr = parMapDFTArray (VS.zipWith (*) bias) convolvedArr
      !maxMag =
        L.sum .
        parMap rdeepseq (VS.sum . VS.map (\x -> (magnitude x) ^ 2)) .
        getDFTArrayVector $
        biasedConvolvedArr
      !normalizedBiasedConvolvedArr =
        parMapDFTArray (VS.map (/ (maxMag :+ 0))) biasedConvolvedArr
  when writeFlag .
    plotDFTArrayPower
      (folderPath </> (printf "Source_%03d.png" numStep))
      cols
      rows $
    convolvedArr
  powerMethod'
    plan
    folderPath
    writeFlag
    coefficients
    harmonicsArray
    bias
    (numStep - 1)
    thetaRHarmonics
    normalizedBiasedConvolvedArr


{-# INLINE computeContour' #-}
computeContour' ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> String -> [[Complex Double]]
  -> DFTArray
  -> IO DFTArray
computeContour' !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !numIteration !suffix thetaRHarmonics !input@(DFTArray rows cols _ _ _) = do
  eigenVector <-
    powerMethod'
      plan
      folderPath
      writeFlag
      coefficients
      harmonicsArray
      bias
      numIteration
      thetaRHarmonics
      input
  source <- convolve' Source plan coefficients harmonicsArray eigenVector
  sink <- convolve' Sink plan coefficients harmonicsArray eigenVector
  completion <- completionField' plan source sink
  plotDFTArrayPower (folderPath </> "Sink.png") rows cols sink
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower_%s.png" suffix))
    rows
    cols
    completion
  let point = (32, -4)
      source' =
        fromUnboxed (Z :. 72 :. cols :. rows) .
        VS.convert .
        VS.concat .
        computeFourierSeriesThetaR thetaRHarmonics .
        getDFTArrayVector $
        source
      sink' =
        fromUnboxed (Z :. 72 :. cols :. rows) .
        VS.convert .
        VS.concat .
        computeFourierSeriesThetaR thetaRHarmonics .
        getDFTArrayVector $
        sink
      completion' =
        fromUnboxed (Z :. 72 :. cols :. rows) .
        VS.convert .
        VS.concat .
        computeFourierSeriesThetaR thetaRHarmonics . getDFTArrayVector $
        completion
  plotThetaDimension folderPath ("FourierSeries_Completion_") point .
    R.map magnitude $
    completion'
  plotThetaDimension
    folderPath
    ("FourierSeries_Source_")
    point .
    R.map magnitude $
    source'
  plotThetaDimension
    folderPath
    ("FourierSeries_Sink_")
    point .
    R.map magnitude $
    sink'
  return completion

powerMethodG' ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> VS.Vector (Complex Double)
  -> Int
  -> DFTArray
  -> IO DFTArray
powerMethodG' _ _ _ _ _ _ _ 0 !arr = return arr
powerMethodG' !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !gaussian !numStep !input@(DFTArray rows cols _ _ _) = do
  printCurrentTime (show numStep)
  convolvedArr <- convolveG' Source plan coefficients harmonicsArray gaussian input
  let !biasedConvolvedArr = parMapDFTArray (VS.zipWith (*) bias) convolvedArr
      !maxMag =
        L.sum .
        parMap rdeepseq (VS.sum . VS.map (\x -> (magnitude x)^2)) .
        getDFTArrayVector $
        biasedConvolvedArr
      !normalizedBiasedConvolvedArr =
        parMapDFTArray (VS.map (/ (maxMag :+ 0))) biasedConvolvedArr
  when writeFlag .
    plotDFTArrayPower
      (folderPath </> (printf "Source_%03d.png" numStep))
      cols
      rows $
    convolvedArr
  powerMethodG'
    plan
    folderPath
    writeFlag
    coefficients
    harmonicsArray
    bias
    gaussian
    (numStep - 1)
    normalizedBiasedConvolvedArr


{-# INLINE computeContourG' #-}
computeContourG' ::
     DFTPlan
  -> FilePath
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> VS.Vector (Complex Double)
  -> Int
  -> String
  -> DFTArray
  -> IO DFTArray
computeContourG' !plan !folderPath !writeFlag !coefficients !harmonicsArray !bias !gaussian !numIteration !suffix !input@(DFTArray rows cols _ _ _) = do
  eigenVector <-
    powerMethodG'
      plan
      folderPath
      writeFlag
      coefficients
      harmonicsArray
      bias
      gaussian
      numIteration
      input
  source <- convolveG' Source plan coefficients harmonicsArray gaussian eigenVector
  sink <- convolveG' Sink plan coefficients harmonicsArray gaussian eigenVector
  completion <- completionField' plan source sink
  plotDFTArrayPower (folderPath </> "Sink.png") rows cols sink
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower_%s.png" suffix))
    rows
    cols
    completion
  return completion

powerMethodSparse ::
     Field
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> [(Double, Double)]
  -> Int
  -> [R.Array U DIM2 (Complex Double)]
  -> [R.Array U DIM2 (Complex Double)]
powerMethodSparse _ _ _ _ _ _ _ 0 !arr = arr
powerMethodSparse !field !coefficients !harmonicsArray !thetaFreqs !rFreqs !cutoff !xs !numStep !input =
  let !convolvedArr =
        convolveSparse
          field
          coefficients
          harmonicsArray
          thetaFreqs
          rFreqs
          cutoff
          xs
          input
      !maxMag =
        sqrt .
        L.sum .
        parMap rdeepseq (VU.sum . VU.map (\x -> (magnitude x) ^ 2) . toUnboxed) $
        convolvedArr
      !normalizedConvolvedArr =
        parMap rseq (computeS . R.map (/ (maxMag :+ 0))) convolvedArr
  in powerMethodSparse
       field
       coefficients
       harmonicsArray
       thetaFreqs
       rFreqs
       cutoff
       xs
       (numStep - 1)
       normalizedConvolvedArr

{-# INLINE computeContourSparse #-}
computeContourSparse ::
     DFTPlan
  -> FilePath
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> [(Double, Double)]
  -> Int
  -> String
  -> [R.Array U DIM2 (Complex Double)]
  -> IO [R.Array U DIM2 (Complex Double)]
  -- -> IO DFTArray
computeContourSparse !plan !folderPath !coefficients !harmonicsArray !thetaFreqs !rFreqs !cutoff !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
      ys =
        powerMethodSparse
          Source
          coefficients
          harmonicsArray
          thetaFreqs
          rFreqs
          cutoff
          xs
          numIteration
          input
      !eigenSourceSparse =
        sparseArrayToDFTArray
          rows
          cols
          (R.toList thetaFreqs)
          (R.toList rFreqs)
          xs
          ys
  harmonicsArrayDFT <-
    fmap (listArray (bounds harmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    harmonicsArray
  printCurrentTime "Source"
  source <-
    STC.Convolution.convolve Source plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SourcePower_Sparse.png") rows cols source
  printCurrentTime "Sink"
  let sink =
        parZipWithDFTArray
          (\vec (rFreq, thetaFreq) -> VS.map (* (cis (thetaFreq * pi))) vec)
          source
          ((,) <$> (R.toList rFreqs) <*> (R.toList thetaFreqs))
  -- sink <- convolve Sink plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SinkPower_Sparse.png") rows cols sink
  printCurrentTime "Completion"
  completion <- completionField plan source sink
  let mm =
        L.maximum . L.map (VS.maximum . VS.map magnitude) . getDFTArrayVector $
        completion
  print mm
  plotDFTArrayPower
    (folderPath </> ("CompletionPower_Sparse.png"))
    rows
    cols .
    parMapDFTArray (VS.map (/ (mm :+ 0))) $
    completion
  return ys

{-# INLINE computeContourSparseG #-}
computeContourSparseG ::
     DFTPlan
  -> FilePath
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> VS.Vector (Complex Double)
  -> [(Double, Double)]
  -> Int
  -> String
  -> [R.Array U DIM2 (Complex Double)]
  -> IO DFTArray
computeContourSparseG !plan !folderPath !coefficients !harmonicsArray !thetaFreqs !rFreqs !cutoff !gaussian !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
  blurredVecF <-
    fmap (L.map (VS.zipWith (*) gaussian)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VU.convert . toUnboxed) . IA.elems $
    harmonicsArray
  blurredHarmonicsArray <-
    (listArray (bounds harmonicsArray) .
     L.map
       (\vec ->
          let arr = fromUnboxed (Z :. cols :. rows) . VS.convert $ vec
          in computeS . R.traverse arr id $ \f idx@(Z :. i :. j) ->
               let r =
                     sqrt . fromIntegral $
                     (i - div cols 2) ^ 2 + (j - div rows 2) ^ 2
               in if r <= 0
                    then 0
                    else f idx)) <$>
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) blurredVecF
  let !eigenSourceSparse =
        sparseArrayToDFTArray
          rows
          cols
          (R.toList thetaFreqs)
          (R.toList rFreqs)
          xs .
        powerMethodSparse
          Source
          coefficients
          blurredHarmonicsArray
          thetaFreqs
          rFreqs
          cutoff
          xs
          numIteration $
        input
  harmonicsArrayDFT <-
    fmap (listArray (bounds blurredHarmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    blurredHarmonicsArray
  printCurrentTime "Source"
  source <-
    STC.Convolution.convolve Source plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SourcePower.png") rows cols source
  printCurrentTime "Sink"
  let sink =
        parZipWithDFTArray
          (\vec (rFreq, thetaFreq) -> VS.map (* (cis (thetaFreq * pi))) vec)
          source
          ((,) <$> (R.toList rFreqs) <*> (R.toList thetaFreqs))
  -- sink <- convolve Sink plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SinkPower.png") rows cols sink
  printCurrentTime "Completion"
  completion <- completionField plan source sink
  let mm =
        L.maximum . L.map (VS.maximum . VS.map magnitude) . getDFTArrayVector $
        completion
  print mm
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower_%s.png" suffix))
    rows
    cols .
    parMapDFTArray (VS.map (/ (mm :+ 0))) $
    completion
  return completion

powerMethodSparse' ::
     Field
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> [(Double, Double)]
  -> Int
  -> [R.Array U DIM1 (Complex Double)]
  -> IO [R.Array U DIM1 (Complex Double)]
powerMethodSparse' _ _ _ _ _ _ _ 0 !arr = return arr
powerMethodSparse' !field !coefficients !harmonicsArray !phiFreqs !rhoFreqs  !cutoff !xs !numStep !input =
  let !convolvedArr =
        convolveSparse'
          field
          coefficients
          harmonicsArray
          phiFreqs
          rhoFreqs
          cutoff
          xs
          input
      !maxMag =
        sqrt . L.sum .
        parMap rdeepseq (VU.sum . VU.map (\x -> (magnitude x) ^ 2) . toUnboxed) $
        convolvedArr
      !normalizedConvolvedArr =
        parMap rseq (computeS . R.map (\x -> x / (maxMag :+ 0))) convolvedArr
  in do
        powerMethodSparse'
          field
          coefficients
          harmonicsArray
          phiFreqs
          rhoFreqs
          cutoff
          xs
          (numStep - 1)
          normalizedConvolvedArr


{-# INLINE computeContourSparse' #-}
computeContourSparse' ::
     DFTPlan
  -> FilePath
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> [[Complex Double]]
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> [(Double, Double)]
  -> Int
  -> String
  -> [R.Array U DIM1 (Complex Double)]
  -> IO [R.Array U DIM1 (Complex Double)]
  -- -> IO DFTArray
computeContourSparse' !plan !folderPath !coefficients !harmonicsArray thetaRHarmonics !phiFreqs !rhoFreqs !cutoff !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
  -- !eigenSourceSparse <-
  --   sparseArrayToDFTArray' rows cols (R.toList phiFreqs) (R.toList rhoFreqs) xs <$>
    -- powerMethodSparse'
    --   Source
    --   coefficients
    --   harmonicsArray
    --   phiFreqs
    --   rhoFreqs
    --   cutoff
    --   xs
    --   numIteration
    --   input
  ys <-
    powerMethodSparse'
      Source
      coefficients
      harmonicsArray
      phiFreqs
      rhoFreqs
      cutoff
      xs
      numIteration
      input
  let eigenSourceSparse =
        sparseArrayToDFTArray'
          rows
          cols
          (R.toList phiFreqs)
          (R.toList rhoFreqs)
          xs
          ys
  print (bounds harmonicsArray)
  print . L.length . IA.elems $ harmonicsArray
  harmonicsArrayDFT <-
    fmap (listArray (bounds harmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    harmonicsArray
  printCurrentTime "Source"
  source <-
    convolve' Source plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SourcePower_Sparse.png") rows cols source
  -- plotDFTArrayThetaR
  --   (folderPath </> "Source.png")
  --   rows
  --   cols
  --   thetaRHarmonics
  --   source
  printCurrentTime "Sink"
  -- sink <- convolve' Sink plan coefficients harmonicsArrayDFT eigenSourceSparse
  let sink =
        parZipWithDFTArray
          (\vec thetaFreq -> VS.map (* (cis (thetaFreq * pi))) vec)
          source
          (R.toList phiFreqs)
  plotDFTArrayPower (folderPath </> "SinkPower_Sparse.png") rows cols sink
  -- plotDFTArrayThetaR (folderPath </> "Sink.png") rows cols thetaRHarmonics sink
  printCurrentTime "Completion"
  completion <- completionField' plan source sink
  let mm =
        L.maximum . L.map (VS.maximum . VS.map magnitude) . getDFTArrayVector $
        completion
  print mm
  plotDFTArrayPower (folderPath </> ("CompletionPower_Sparse.png")) rows cols -- .
    -- parMapDFTArray (VS.map (/ (mm :+ 0))) $
    completion
  -- let point = (31, -4)
  --     source' =
  --       fromUnboxed (Z :. 72 :. cols :. rows) .
  --       VS.convert .
  --       VS.concat .
  --       computeFourierSeriesThetaR thetaRHarmonics .
  --       getDFTArrayVector $
  --       source
  --     sink' =
  --       fromUnboxed (Z :. 72 :. cols :. rows) .
  --       VS.convert .
  --       VS.concat .
  --       computeFourierSeriesThetaR thetaRHarmonics .
  --       getDFTArrayVector $
  --       sink
  --     completion' =
  --       fromUnboxed (Z :. 72 :. cols :. rows) .
  --       VS.convert .
  --       VS.concat .
  --       computeFourierSeriesThetaR thetaRHarmonics .
  --       getDFTArrayVector $
  --       completion
  -- plotThetaDimension
  --   folderPath
  --   ("FourierSeries_Completion_")
  --   point .
  --   R.map magnitude $
  --   completion'
  -- plotThetaDimension
  --   folderPath
  --   ("FourierSeries_Source_")
  --   point .
  --   R.map magnitude $
  --   source'
  -- plotThetaDimension
  --   folderPath
  --   ("FourierSeries_Sink_")
  --   point .
  --   R.map magnitude $
  --   sink'
  -- return completion
  return ys

{-# INLINE computeContourSparseG' #-}
computeContourSparseG' ::
     DFTPlan
  -> FilePath
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> [[Complex Double]]
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> VS.Vector (Complex Double)
  -> Double
  -> [(Double, Double)]
  -> Int
  -> String
  -> [R.Array U DIM1 (Complex Double)]
  -> IO DFTArray
computeContourSparseG' !plan !folderPath !coefficients !harmonicsArray !thetaRHarmonics !phiFreqs !rhoFreqs !cutoff !gaussian !std !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
  blurredVecF <-
    fmap (L.map (VS.zipWith (*) gaussian)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VU.convert . toUnboxed) . IA.elems $
    harmonicsArray
  blurredHarmonicsArray <-
    (listArray (bounds harmonicsArray) .
     L.map
       (\vec ->
          let arr = fromUnboxed (Z :. cols :. rows) . VS.convert $ vec
          in computeS . R.traverse arr id $ \f idx@(Z :. i :. j) ->
               let r =
                     sqrt . fromIntegral $
                     (i - div cols 2) ^ 2 + (j - div rows 2) ^ 2
               in if r <= 0
                    then 0
                    else f idx)) <$>
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) blurredVecF
  !eigenSourceSparse <-
    sparseArrayToDFTArray'
      rows
      cols
      (R.toList phiFreqs)
      (R.toList rhoFreqs)
      xs <$>
    powerMethodSparse'
      Source
      coefficients
      blurredHarmonicsArray
      phiFreqs
      rhoFreqs
      cutoff
      xs
      numIteration
      input
  harmonicsArrayDFT <-
    fmap (listArray (bounds blurredHarmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    blurredHarmonicsArray
  printCurrentTime "Source"
  source <-
    convolve' Source plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SourcePower.png") rows cols source
  plotDFTArrayThetaR
    (folderPath </> "Source.png")
    rows
    cols
    thetaRHarmonics
    source
  printCurrentTime "Sink"
  sink <- convolve' Sink plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SinkPower.png") rows cols sink
  plotDFTArrayThetaR (folderPath </> "Sink.png") rows cols thetaRHarmonics sink
  printCurrentTime "Completion"
  completion <- completionField' plan source sink
  let mm =
        L.maximum . L.map (VS.maximum . VS.map magnitude) . getDFTArrayVector $
        completion
  print mm
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower_%s.png" suffix))
    rows
    cols .
    parMapDFTArray (VS.map (/ (mm :+ 0))) $
    completion
  return completion


-- # INLINE computeContourSparse''' #
computeContourSparse''' ::
     DFTPlan
  -> FilePath
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> [[Complex Double]]
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double -> VS.Vector (Complex Double)
  -> [(Double, Double)]
  -> Int
  -> String
  -> [R.Array U DIM1 (Complex Double)]
  -> IO DFTArray
computeContourSparse''' !plan !folderPath !coefficients !harmonicsArray !thetaRHarmonics !phiFreqs !rhoFreqs !cutoff !gaussian !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
  blurredVecF <-
    fmap (L.map (VS.zipWith (*) gaussian)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VU.convert . toUnboxed) . IA.elems $
    harmonicsArray
  blurredHarmonicsArray <-
    (listArray (bounds harmonicsArray) .
     L.map
       (\vec ->
          let arr = fromUnboxed (Z :. cols :. rows) . VS.convert $ vec
          in computeS . R.traverse arr id $ \f idx@(Z :. i :. j) ->
               if (i == div cols 2) && (j == div rows 2)
                 then 0
                 else f idx)) <$>
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) blurredVecF
  !eigenSourceSparse <-
    sparseArrayToDFTArray'''
      rows
      cols
      (R.toList phiFreqs)
      (R.toList rhoFreqs)
      xs <$>
    powerMethodSparse'
      Source
      coefficients
      blurredHarmonicsArray
      phiFreqs
      rhoFreqs
      cutoff
      xs
      numIteration
      input
  harmonicsArrayDFT <-
    fmap (listArray (bounds blurredHarmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    blurredHarmonicsArray
  xs <-
    M.mapM
      (\(i, (p, q)) -> do
         source <- convolve' Source plan coefficients harmonicsArrayDFT q
         plotDFTArrayPower
           (folderPath </> (printf "SourcePower_%d.png" i))
           rows
           cols
           source
         sink <- convolve' Sink plan coefficients harmonicsArrayDFT p
         plotDFTArrayPower
           (folderPath </> (printf "SinkPower_%d.png" i))
           rows
           cols
           sink
         completion <- completionField' plan source sink
         plotDFTArrayPower
           (folderPath </> (printf "CompletionPower_%d.png" i))
           rows
           cols
           completion
         when
           (i == 0)
           (do let point = (31, -4)
                   source' =
                     fromUnboxed (Z :. 72 :. cols :. rows) .
                     VS.convert .
                     VS.concat .
                     computeFourierSeriesThetaR thetaRHarmonics .
                     getDFTArrayVector $
                     source
                   sink' =
                     fromUnboxed (Z :. 72 :. cols :. rows) .
                     VS.convert .
                     VS.concat .
                     computeFourierSeriesThetaR thetaRHarmonics .
                     getDFTArrayVector $
                     sink
                   completion' =
                     fromUnboxed (Z :. 72 :. cols :. rows) .
                     VS.convert .
                     VS.concat .
                     computeFourierSeriesThetaR thetaRHarmonics .
                     getDFTArrayVector $
                     completion
               plotThetaDimension
                 folderPath
                 ("FourierSeries_Completion_")
                 point .
                 R.map magnitude $
                 completion'
               plotThetaDimension
                 folderPath
                 ("FourierSeries_Source_")
                 point .
                 R.map magnitude $
                 source'
               plotThetaDimension
                 folderPath
                 ("FourierSeries_Sink_")
                 point .
                 R.map magnitude $
                 sink')
         return completion) .
    L.zip [0 :: Int ..] $
    eigenSourceSparse
  let (DFTArray rows cols thetaFreqs rFreqs _) = L.head xs
      completion =
        DFTArray rows cols thetaFreqs rFreqs .
        L.foldl1' (L.zipWith (VS.zipWith (+))) .
        L.map
          (\arr ->
             let ys = getDFTArrayVector arr
                 s = L.maximum . L.map (VS.maximum . VS.map magnitude) $ ys
             in L.map (VS.map (/ (s :+ 0))) ys) $
        xs
      mm =
        L.maximum . L.map (VS.maximum . VS.map magnitude) . getDFTArrayVector $
        completion
  print mm
  plotDFTArrayPower
    (folderPath </> (printf "CompletionPower_%s.png" suffix))
    rows
    cols .
    parMapDFTArray (VS.map (/ (mm :+ 0))) $
    completion
  return completion

powerMethodPinwheelBasis' ::
     DFTPlan
  -> FilePath
  -> [PTX]
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> Int
  -> DFTArray
  -> IO DFTArray
powerMethodPinwheelBasis' _ _ _ _ _ _ _ _ _ _ _ _ 0 arr = return arr
powerMethodPinwheelBasis' plan folderPath ptxs writeFlag coefficients harmonicsArray dftBias numPoints numR2Freqs delta periodR2 numBatch numStep input@(DFTArray rows cols _ _ _) = do
  printCurrentTime (show numStep)
  let convolvedArr = convoluvePinhweelBasis' coefficients harmonicsArray input
  biasedConvolvedArr <-
    multiplyPinwheelBasisBatch1 plan numBatch dftBias convolvedArr
  let s =
        L.maximum .
        parMap rdeepseq (VS.maximum . VS.map (\x -> (magnitude x) )) .
        getDFTArrayVector $
        biasedConvolvedArr
      normalizedBiasedConvolvedArr =
        parMapDFTArray (VS.map (/ (s :+ 0))) biasedConvolvedArr
      -- envelope =
      --   VU.convert .
      --   toUnboxed . computeS . fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
      --     let r =
      --           sqrt . fromIntegral $
      --           (i - div numR2Freqs 2) ^ 2 + (j - div numR2Freqs 2) ^ 2
      --     in if r == 0
      --          then 0
      --          else r ** (-0.75) :+ 0
      -- normalizedBiasedConvolvedArr =
      --   parMapDFTArray (VS.zipWith (*) envelope) normalizedBiasedConvolvedArr'
  when
    writeFlag
    (do let eigenMat =
              A.transpose . A.use . fromDFTArray . getDFTArrayVector $
              convolvedArr
        eigenR2 <-
          computeFourierSeriesR2StreamAcc
            ptxs
            numR2Freqs
            numPoints
            (L.length . getDFTArrayVector $ convolvedArr)
            periodR2
            delta
            numBatch
            eigenMat
        plotImageRepa (folderPath </> (printf "Source_%03d.png" numStep)) .
          ImageRepa 8 .
          fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
          -- VU.map sqrt .
          toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
          eigenR2
        -- let biasMat =
        --       A.transpose . A.use . fromDFTArray . getDFTArrayVector $
        --       normalizedBiasedConvolvedArr
        -- biasR2 <-
        --   computeFourierSeriesR2StreamAcc
        --     ptxs
        --     numR2Freqs
        --     numPoints
        --     (L.length . getDFTArrayVector $ convolvedArr)
        --     periodR2
        --     delta
        --     1
        --     biasMat
        -- plotImageRepa (folderPath </> (printf "Biased_%03d.png" numStep)) .
        --   ImageRepa 8 .
        --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
        --   VU.map sqrt .
        --   toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
        --   biasR2
    )
  powerMethodPinwheelBasis'
    plan
    folderPath
    ptxs
    writeFlag
    coefficients
    harmonicsArray
    dftBias
    numPoints
    numR2Freqs
    delta
    periodR2
    numBatch
    (numStep - 1)
    normalizedBiasedConvolvedArr


computeContourPinwheelBasis' ::
     DFTPlan
  -> FilePath
  -> [PTX]
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> Int
  -> DFTArray
  -> IO (R.Array U DIM4 (Complex Double))
computeContourPinwheelBasis' plan folderPath ptxs writeFlag coefficients harmonicsArray dftBias numStep numBatch numPoints numR2Freqs delta periodR2 thetaFreqs rFreqs input = do
  eigenSource' <-
    powerMethodPinwheelBasis'
      plan
      folderPath
      ptxs
      writeFlag
      coefficients
      harmonicsArray
      dftBias
      numPoints
      numR2Freqs
      delta
      periodR2
      numBatch
      numStep
      input
  let eigenSource =
        convoluvePinhweelBasis' coefficients harmonicsArray eigenSource'
      numLogPolarFreqs = thetaFreqs * rFreqs
      eigenSourceMat =
        A.transpose . A.use . fromDFTArray . getDFTArrayVector $ eigenSource
  eigenSourceR2' <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freqs
      numPoints
      numLogPolarFreqs
      periodR2
      delta
      numBatch
      eigenSourceMat
  let eigenSourceR2 =
        R.reshape (Z :. rFreqs :. thetaFreqs :. numPoints :. numPoints) $
        eigenSourceR2'
      eigenSinkR2 =
        timeReversalRepa
          (L.map fromIntegral . getListFromNumber $ thetaFreqs)
          eigenSourceR2
  completionR2 <- completionFieldRepa plan eigenSourceR2 eigenSinkR2
  plotImageRepa (folderPath </> "Source.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    -- VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    eigenSourceR2
  plotImageRepa (folderPath </> "Sink.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    -- VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    eigenSinkR2
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    -- VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    completionR2
  return completionR2
  
powerMethodPinwheelBasis ::
     DFTPlan
  -> FilePath
  -> [PTX]
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> Int
  -> DFTArray
  -> IO DFTArray
powerMethodPinwheelBasis _ _ _ _ _ _ _ _ _ _ _ _ 0 arr = return arr
powerMethodPinwheelBasis plan folderPath ptxs writeFlag coefficients harmonicsArray dftBias numPoints numR2Freqs delta periodR2 numBatch numStep input@(DFTArray rows cols _ _ _) = do
  printCurrentTime (show numStep)
  let convolvedArr = convoluvePinhweelBasis coefficients harmonicsArray input
  biasedConvolvedArr <-
    multiplyPinwheelBasisBatch1 plan numBatch dftBias convolvedArr
  let s =
        L.maximum .
        parMap rdeepseq (VS.maximum . VS.map (\x -> (magnitude x) )) .
        getDFTArrayVector $
        biasedConvolvedArr
      normalizedBiasedConvolvedArr =
        parMapDFTArray (VS.map (/ (s :+ 0))) biasedConvolvedArr
      -- envelope =
      --   VU.convert .
      --   toUnboxed . computeS . fromFunction (Z :. numR2Freqs :. numR2Freqs) $ \(Z :. i :. j) ->
      --     let r =
      --           sqrt . fromIntegral $
      --           (i - div numR2Freqs 2) ^ 2 + (j - div numR2Freqs 2) ^ 2
      --     in if r == 0
      --          then 0
      --          else r ** (-0.75) :+ 0
      -- normalizedBiasedConvolvedArr =
      --   parMapDFTArray (VS.zipWith (*) envelope) normalizedBiasedConvolvedArr'
  when
    writeFlag
    (do let eigenMat =
              A.transpose . A.use . fromDFTArray . getDFTArrayVector $
              convolvedArr
        eigenR2 <-
          computeFourierSeriesR2StreamAcc
            ptxs
            numR2Freqs
            numPoints
            (L.length . getDFTArrayVector $ convolvedArr)
            periodR2
            delta
            numBatch
            eigenMat
        plotImageRepa (folderPath </> (printf "Source_%03d.png" numStep)) .
          ImageRepa 8 .
          fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
          VU.map sqrt .
          toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
          eigenR2
        let biasMat =
              A.transpose . A.use . fromDFTArray . getDFTArrayVector $
              normalizedBiasedConvolvedArr
        biasR2 <-
          computeFourierSeriesR2StreamAcc
            ptxs
            numR2Freqs
            numPoints
            (L.length . getDFTArrayVector $ convolvedArr)
            periodR2
            delta
            numBatch
            biasMat
        plotImageRepa (folderPath </> (printf "Biased_%03d.png" numStep)) .
          ImageRepa 8 .
          fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
          VU.map sqrt .
          toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
          biasR2
    )
  powerMethodPinwheelBasis
    plan
    folderPath
    ptxs
    writeFlag
    coefficients
    harmonicsArray
    dftBias
    numPoints
    numR2Freqs
    delta
    periodR2
    numBatch
    (numStep - 1)
    normalizedBiasedConvolvedArr


computeContourPinwheelBasis ::
     DFTPlan
  -> FilePath
  -> [PTX]
  -> Bool
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> Int
  -> Int
  -> DFTArray
  -> IO (R.Array U DIM4 (Complex Double))
computeContourPinwheelBasis plan folderPath ptxs writeFlag coefficients harmonicsArray dftBias numStep numBatch numPoints numR2Freqs delta periodR2 thetaFreqs rFreqs input = do
  eigenSource' <-
    powerMethodPinwheelBasis
      plan
      folderPath
      ptxs
      writeFlag
      coefficients
      harmonicsArray
      dftBias
      numPoints
      numR2Freqs
      delta
      periodR2
      numBatch
      numStep
      input
  let eigenSource =
        convoluvePinhweelBasis coefficients harmonicsArray eigenSource'
      numLogPolarFreqs = thetaFreqs * rFreqs
      eigenSourceMat =
        A.transpose . A.use . fromDFTArray . getDFTArrayVector $ eigenSource
  eigenSourceR2' <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freqs
      numPoints
      numLogPolarFreqs
      periodR2
      delta
      numBatch
      eigenSourceMat
  let eigenSourceR2 =
        R.reshape (Z :. rFreqs :. thetaFreqs :. numPoints :. numPoints) $
        eigenSourceR2'
      eigenSinkR2 =
        timeReversalRepa
          (L.map fromIntegral . getListFromNumber $ thetaFreqs)
          eigenSourceR2
  completionR2 <- completionFieldRepa plan eigenSourceR2 eigenSinkR2
  plotImageRepa (folderPath </> "Source.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    eigenSourceR2
  plotImageRepa (folderPath </> "Sink.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    eigenSinkR2
  plotImageRepa (folderPath </> "Completion.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    completionR2
  return completionR2



powerMethodFourierPinwheel ::
     DFTPlan
  -> FilePath
  -> Bool
  -> FPData (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Double
  -> Double
  -> Int -> VS.Vector (Complex Double)
  -> Int 
  -> FPArray (VS.Vector (Complex Double))
  -> IO (FPArray (VS.Vector (Complex Double)))
powerMethodFourierPinwheel _ _ _ _ _ _ _ _ _ _ 0 arr = return arr
powerMethodFourierPinwheel plan folderPath writeFlag harmonicsArray dftBias numPoints delta periodR2 numBatch gaussianEnvelope numStep input = do
  printCurrentTime (show numStep)
  convolvedArr <- FP.convolve harmonicsArray input
  biasedConvolvedArr <- multiplyBias4D plan dftBias convolvedArr
  let s =
        L.sum .
        parMap rdeepseq (VS.sum . VS.map (\x -> (magnitude x) ^ 2)) . getFPArray $
        biasedConvolvedArr
      normalizedBiasedConvolvedArr =
        parMapFPArray (VS.map (/ (s :+ 0))) biasedConvolvedArr
      -- normalizedBiasedConvolvedArr = parMapFPArray (VS.zipWith (*) gaussianEnvelope) normalizedBiasedConvolvedArr'
  when
    (writeFlag && mod numStep 2 == 0)
    (do _ <-
          plotFPArray
            plan
            (folderPath </> (printf "Source_%03d.png" numStep))
            convolvedArr
        -- _ <- plotFPArray
        --        plan
        --        (folderPath </> (printf "Bias_%03d.png" numStep))
        --        normalizedBiasedConvolvedArr
        -- plotFPArrayFreqency
        --   (folderPath </> (printf "SourceFreq_%03d.png" numStep))
        --   convolvedArr
        -- plotFPArrayFreqency
        --   (folderPath </> (printf "BiasFreq1_%03d.png" numStep))
        --   normalizedBiasedConvolvedArr' 
        -- plotFPArrayFreqency
        --   (folderPath </> (printf "BiasFreq2_%03d.png" numStep))
        --   normalizedBiasedConvolvedArr 
        return ())
  powerMethodFourierPinwheel
    plan
    folderPath
    writeFlag
    harmonicsArray
    dftBias
    numPoints
    delta
    periodR2
    numBatch
    gaussianEnvelope
    (numStep - 1)
    normalizedBiasedConvolvedArr


computeContourFourierPinwheel ::
     DFTPlan
  -> FilePath
  -> Bool
  -> FPData (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double -> Double
  -> FPArray (VS.Vector (Complex Double))
  -> String
  -> [Int]
  -> VS.Vector (Complex Double)
  -> IO (R.Array U DIM4 (Complex Double))
computeContourFourierPinwheel plan folderPath writeFlag harmonicsArray dftBias numStep numBatch numPoints delta periodR2 periodEnv input suffix deviceIDs gaussianEnvelope = do
  eigenSource' <-
    powerMethodFourierPinwheel
      plan
      folderPath
      writeFlag
      harmonicsArray
      dftBias
      numPoints
      delta
      periodR2
      numBatch
      gaussianEnvelope
      numStep
      input
  eigenSource <- FP.convolve harmonicsArray eigenSource'
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  -- eigenSourceR2 <-
  --   plotFPArray plan (folderPath </> printf "Source_%s.png" suffix) eigenSource
  eigenSourceR2 <-
    plotFPArrayAcc
      ptxs
      (folderPath </> printf "Source_%s.png" suffix)
      numPoints
      delta
      periodR2
      numBatch
      eigenSource
  let numRFreq = getFPArrayNumRFreq eigenSource
      numThetaFreq = getFPArrayNumThetaFreq eigenSource
  eigenSinkR2 <-
    computeUnboxedP $
    timeReversalRepa
      (L.map fromIntegral . getListFromNumber $ numThetaFreq)
      eigenSourceR2
  completionR2 <- completionFieldRepa plan eigenSourceR2 eigenSinkR2
  -- plotImageRepa (folderPath </> printf "Source_%s.png" suffix) .
  --   ImageRepa 8 .
  --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
  --   -- VU.map (\x -> x^2) .
  --   toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
  --   eigenSourceR2
  (sumP . R.map (\x -> magnitude x ** 2) . rotate4D2 $ eigenSinkR2) >>= sumP >>=
    plotImageRepa (folderPath </> printf "Sink_%s.png" suffix) .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All)
  completionR2' <-
    (sumP . R.map (\x -> magnitude x ** 2) . rotate4D2 $ completionR2) >>= sumP
  plotImageRepa (folderPath </> printf "Completion_%s.png" suffix) .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    completionR2'
  plotImageRepa (folderPath </> printf "Completion1_%s.png" suffix) .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) . R.map sqrt $
    completionR2'
  return completionR2


powerMethodFourierPinwheelDiscrete ::
     DFTPlan
  -> FilePath
  -> Bool
  -> FPData (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Double
  -> Double
  -> Int
  -> (Int, Int)
  -> Int
  -> FPArray (VS.Vector (Complex Double))
  -> IO (FPArray (VS.Vector (Complex Double)))
powerMethodFourierPinwheelDiscrete  _ _ _ _ _ _ _ _ _ idx 0 arr = return arr
powerMethodFourierPinwheelDiscrete plan folderPath writeFlag harmonicsArray dftBias numPoints delta periodR2 numBatch idx numStep input' = do
  printCurrentTime (show numStep)
  -- input <- multiplyRFunction plan (periodR2^2 / 2) input'
  convolvedArr <- FP.convolve harmonicsArray input'
  biasedConvolvedArr <- multiplyBiasDiscrete plan dftBias convolvedArr
  let s =
        L.sum .
        parMap rdeepseq (VS.sum . VS.map (\x -> (magnitude x) ^ 2)) . getFPArray $
        biasedConvolvedArr
      normalizedBiasedConvolvedArr =
        parMapFPArray (VS.map (/ (s :+ 0))) biasedConvolvedArr
  when
    (writeFlag && mod numStep 2 == 0
    )
    (do arrR2 <-
          plotFPArray
            plan
            (folderPath </> (printf "Source_%03d.png" numStep))
            convolvedArr
        plotRThetaDist
          (folderPath </> (printf "ThetaDist_%03d.png" numStep))
          (folderPath </> (printf "RDist_%03d.png" numStep))
          (getFPArrayNumXFreq convolvedArr)
          72
          100
          periodR2
          idx
          arrR2
        _ <- plotFPArray
               plan
               (folderPath </> (printf "Bias_%03d.png" numStep))
               normalizedBiasedConvolvedArr
        return ()
        -- let biasMat = A.transpose . toMatrixAcc $ normalizedBiasedConvolvedArr
        -- biasR2 <-
        --   computeFourierSeriesR2StreamAcc
        --     ptxs
        --     (getFPArrayNumXFreq convolvedArr)
        --     numPoints
        --     (getFPArrayNumRFreq convolvedArr *
        --      getFPArrayNumThetaFreq convolvedArr)
        --     periodR2
        --     delta
        --     numBatch
        --     biasMat
        -- plotImageRepa (folderPath </> (printf "Biased_%03d.png" numStep)) .
        --   ImageRepa 8 .
        --   fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
        --   VU.map sqrt .
        --   toUnboxed . sumS . R.map (\x -> (magnitude x) ** 2) . rotate3D $
        --   biasR2
     )
  powerMethodFourierPinwheelDiscrete
    plan
    folderPath
    writeFlag
    harmonicsArray
    dftBias
    numPoints
    delta
    periodR2
    numBatch
    idx
    (numStep - 1)
    normalizedBiasedConvolvedArr


computeContourFourierPinwheelDiscrete ::
     DFTPlan
  -> FilePath
  -> Bool
  -> FPData (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> Int
  -> Int
  -> Int
  -> Double
  -> Double
  -> FPArray (VS.Vector (Complex Double))
  -> String
  -> (Int, Int)
  -> IO (R.Array U DIM4 (Complex Double))
computeContourFourierPinwheelDiscrete plan folderPath writeFlag harmonicsArray dftBias numStep numBatch numPoints delta periodR2 input suffix idx = do
  eigenSource' <-
    powerMethodFourierPinwheelDiscrete
      plan
      folderPath
      writeFlag
      harmonicsArray
      dftBias
      numPoints
      delta
      periodR2
      numBatch
      idx
      numStep
      input
  eigenSource <- FP.convolve harmonicsArray eigenSource'
  eigenSourceR2 <- plotFPArray plan (folderPath </> printf "Source_%s.png" suffix)  eigenSource
  let numThetaFreq = getFPArrayNumThetaFreq eigenSource
      eigenSinkR2 =
        timeReversalRepa
          (L.map fromIntegral . getListFromNumber $ numThetaFreq)
          eigenSourceR2
  completionR2 <- completionFieldRepa plan eigenSourceR2 eigenSinkR2
  plotImageRepa (folderPath </> printf "Sink_%s.png" suffix) .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    eigenSinkR2
  plotImageRepa (folderPath </> printf "Completion_%s.png" suffix) .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numPoints :. numPoints) .
    VU.map sqrt .
    toUnboxed . sumS . sumS . R.map (\x -> (magnitude x) ** 2) . rotate4D2 $
    completionR2
  return completionR2
