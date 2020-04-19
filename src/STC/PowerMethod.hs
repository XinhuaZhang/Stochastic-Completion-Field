{-# LANGUAGE BangPatterns #-}
module STC.PowerMethod where

import           Control.Monad        as M
import           Data.Array.IArray    as IA
import           Data.Array.Repa      as R
import           Data.Complex
import           Data.List            as L
import           Data.Vector.Storable as VS
import           Data.Vector.Unboxed  as VU
import           DFT.Plan
import           Filter.Utils
import           Image.IO
import           STC.CompletionField
import           STC.Convolution
import           STC.DFTArray
import           STC.Plan
import           System.FilePath      ((</>))
import           Text.Printf
import           Utils.Parallel
import           Utils.Time
import Graphics.Gnuplot.Simple
import           FokkerPlanck.FourierSeries

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
    convolve Source plan coefficients harmonicsArray input
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
  source <- convolve Source plan coefficients harmonicsArray eigenVector
  sink <- convolve Sink plan coefficients harmonicsArray eigenVector
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
  -> [(Int, Int)]
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
  -> [(Int, Int)]
  -> Int
  -> String
  -> [R.Array U DIM2 (Complex Double)]
  -> IO DFTArray
computeContourSparse !plan !folderPath !coefficients !harmonicsArray !thetaFreqs !rFreqs !cutoff !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
      !eigenSourceSparse =
        sparseArrayToDFTArray
          rows
          cols
          (R.toList thetaFreqs)
          (R.toList rFreqs)
          xs .
        powerMethodSparse
          Source
          coefficients
          harmonicsArray
          thetaFreqs
          rFreqs
          cutoff
          xs
          numIteration $
        input
  harmonicsArrayDFT <-
    fmap (listArray (bounds harmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    harmonicsArray
  printCurrentTime "Source"
  source <-
    convolve Source plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SourcePower.png") rows cols source
  printCurrentTime "Sink"
  sink <- convolve Sink plan coefficients harmonicsArrayDFT eigenSourceSparse
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
  -> [(Int, Int)]
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
    convolve Source plan coefficients harmonicsArrayDFT eigenSourceSparse
  plotDFTArrayPower (folderPath </> "SourcePower.png") rows cols source
  printCurrentTime "Sink"
  sink <- convolve Sink plan coefficients harmonicsArrayDFT eigenSourceSparse
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
  -> IO DFTArray
computeContourSparse' !plan !folderPath !coefficients !harmonicsArray !thetaRHarmonics !phiFreqs !rhoFreqs !cutoff !xs !numIteration !suffix !input = do
  printCurrentTime $ printf "Power Method %d iterations" numIteration
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
  !eigenSourceSparse <-
    sparseArrayToDFTArray' rows cols (R.toList phiFreqs) (R.toList rhoFreqs) xs <$>
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
  let point = (31, -4)
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
    sink'
  return completion
  
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
    sparseArrayToDFTArrayG'
      2
      std
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
    fmap (listArray (bounds harmonicsArray)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) .
    L.map (VS.convert . toUnboxed . computeUnboxedS . makeFilter2D) . IA.elems $
    harmonicsArray
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


-- {-# INLINE computeContourSparse''' #-}
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
