{-# LANGUAGE BangPatterns #-}
module FourierSeriesOfPinwheels where

import           Control.Monad.Parallel          as MP
import           Data.Array.IArray               as IA
import           Data.Array.Repa                 as R
import           Data.Array.Repa.IO.Binary
import           Data.Array.Repa.Repr.ForeignPtr
import           Data.Complex
import           Data.List                       as L
import           Data.Vector.Storable            as VS
import           Data.Vector.Unboxed             as VU
import           FokkerPlanck.GreensFunction
import           FokkerPlanck.Histogram
import           Image.IO
import           Pinwheel.Base
import           Pinwheel.FourierSeries2D
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Parallel
import           Utils.Time

main = do
  args@(gpuIDsStr:lenStr:deltaStr:inverseDeltaStr:numR2FreqsStr:periodStr:phiFreqsStr:rhoFreqsStr:thetaFreqsStr:scaleFreqsStr:coefFilePath:stdStr:batchSizePolarFreqsStr:batchSizeR2FreqsStr:batchSizeR2Str:numThreadStr:_) <-
    getArgs
  let gpuIDs = read gpuIDsStr :: [Int]
      len = read lenStr :: Double
      delta = read deltaStr :: Double
      inverseDelta = read inverseDeltaStr :: Double
      numR2Freqs = read numR2FreqsStr :: Int
      period = read periodStr :: Double
      phiFreq = read phiFreqsStr :: Int
      phiFreqs = L.map fromIntegral [-phiFreq .. phiFreq]
      rhoFreq = read rhoFreqsStr :: Int
      rhoFreqs = L.map fromIntegral [-rhoFreq .. rhoFreq]
      thetaFreq = read thetaFreqsStr :: Int
      thetaFreqs = L.map fromIntegral [-thetaFreq .. thetaFreq]
      scaleFreq = read scaleFreqsStr :: Int
      scaleFreqs = L.map fromIntegral [-scaleFreq .. scaleFreq]
      numThread = read numThreadStr :: Int
      folderPath = "output/test/FourierSeriesOfPinwheels"
      std = read stdStr :: Double
      batchSizePolarFreqs = read batchSizePolarFreqsStr :: Int
      batchSizeR2Freqs = read batchSizeR2FreqsStr :: Int
      batchSizeR2 = read batchSizeR2Str :: Int
  removePathForcibly (folderPath </> "Recons")
  createDirectoryIfMissing True folderPath
  let !numAngularFreqs = getNumFreqs phiFreq thetaFreq
      !angularFreqsCenter = div numAngularFreqs 2
      !numRadialFreqs = getNumFreqs rhoFreq scaleFreq
      !numPoints' = round $ len / delta :: Int
      !numPoints =
        if odd numPoints'
          then numPoints'
          else numPoints' + 1
  coefficients <-
    computeFourierCoefficients
      gpuIDs
      numR2Freqs
      numPoints
      period
      delta
      batchSizeR2Freqs
      batchSizePolarFreqs
      numAngularFreqs
      numRadialFreqs
  writeArrayToStorableFile coefFilePath coefficients  
  coefficientsArr' <-
    readArrayFromStorableFile
      coefFilePath
      (Z :. numRadialFreqs :. numAngularFreqs :. numR2Freqs :. numR2Freqs) :: IO (R.Array F DIM4 (Complex Double))
  coefficientsArr <- applyGaussian std coefficientsArr'
  let coefficientsVecs =
        L.map
          (\(!rf, !af) ->
             VU.convert . toUnboxed . computeS . R.slice coefficientsArr $
             (Z :. rf :. af :. All :. All))
          [ (rf, af)
          | rf <- [0 .. numRadialFreqs - 1]
          , af <- [0 .. numAngularFreqs - 1]
          ]
  createDirectoryIfMissing True (folderPath </> "Coefficients")
  MP.mapM_
    (\(rf, af) ->
       plotImageRepaComplex
         (folderPath </> "Coefficients" </>
          printf
            "Coef(%d,%d).png"
            (rf - div numRadialFreqs 2)
            (af - div numAngularFreqs 2)) .
       ImageRepa 8 .
       computeS .
       extend (Z :. (1 :: Int) :. All :. All) . R.slice coefficients $
       (Z :. rf :. af :. All :. All))
    [ (rf, af)
    | rf <- [0 .. numRadialFreqs - 1]
    , af <- [0 .. numAngularFreqs - 1]
    ]
  recons <-
    computeFourierSeries
      gpuIDs
      numR2Freqs
      (round len)
      period
      inverseDelta
      batchSizePolarFreqs
      batchSizeR2
      numAngularFreqs
      numRadialFreqs
      coefficientsVecs
  createDirectoryIfMissing True (folderPath </> "Recons")
  MP.mapM_
    (\(!rf, !af) ->
       plotImageRepaComplex
         (folderPath </> "Recons" </>
          printf
            "Recon(%d,%d).png"
            (rf - div numRadialFreqs 2)
            (af - div numAngularFreqs 2)) .
       ImageRepa 8 .
       computeS . extend (Z :. (1 :: Int) :. All :. All) . R.slice recons $
       (Z :. rf :. af :. All :. All)) $
    [ (rf, af)
    | rf <- [0 .. numRadialFreqs - 1]
    , af <- [0 .. numAngularFreqs - 1]
    ]
