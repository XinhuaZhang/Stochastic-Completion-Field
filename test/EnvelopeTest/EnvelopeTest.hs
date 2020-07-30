module EnvelopeTest where

import           Control.Monad                  as M
import           Control.Monad.Parallel         as MP
import qualified Data.Array.Accelerate          as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                as R
import           Data.Binary
import           Data.Complex
import           Data.List                      as L
import           Data.Vector.Storable           as VS
import           Data.Vector.Unboxed            as VU
import           Filter.Utils
import           FokkerPlanck.FourierSeries
import           Foreign.CUDA.Driver            as CUDA
import           FourierMethod.BlockMatrixAcc
import           FourierMethod.FourierSeries2D
import           Image.IO
import           Pinwheel.FourierSeries2D
import           STC
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Array
import           Utils.List
import           Utils.Time

main = do
  args@(deviceIDsStr:numPointsStr:deltaStr:numPointsReconStr:deltaReconStr:numR2FreqStr:periodR2Str:phiFreqsStr:rhoFreqsStr:stdR2Str:stdStr:numBatchR2Str:numBatchR2FreqsStr:sStr:numThreadStr:_) <-
    getArgs
  let deviceIDs = read deviceIDsStr :: [Int]
      numPoints = read numPointsStr :: Int
      delta = read deltaStr :: Double
      numPointsRecon = read numPointsReconStr :: Int
      deltaRecon = read deltaReconStr :: Double
      numR2Freq = read numR2FreqStr :: Int
      periodR2 = read periodR2Str :: Double
      phiFreq = read phiFreqsStr :: Int
      rhoFreq = read rhoFreqsStr :: Int
      numThread = read numThreadStr :: Int
      stdR2 = read stdR2Str :: Double
      std = read stdStr :: Double
      numBatchR2 = read numBatchR2Str :: Int
      numBatchR2Freqs = read numBatchR2FreqsStr :: Int
      s = read sStr :: Double
      folderPath = "output/test/EnvelopeTest"
      func a =
        R.fromFunction (Z :. (1 :: Int) :. numPointsRecon :. numPointsRecon) $ \(Z :. _ :. i :. j) ->
          let r =
                sqrt
                  (fromIntegral $
                   (i - div numPointsRecon 2) ^ 2 +
                   (j - div numPointsRecon 2) ^ 2) :: Double
          in if r == 0
               then 0
               else r ** a :+ 0
  createDirectoryIfMissing True folderPath
  plotImageRepaComplex (folderPath </> "Reference(-0.5).png") .
    ImageRepa 8 . computeS . func $
    (-0.5)
  plotImageRepaComplex (folderPath </> "Reference(-1).png") .
    ImageRepa 8 . computeS . func $
    (-1)
  plotImageRepaComplex (folderPath </> "Reference(-1.5).png") .
    ImageRepa 8 . computeS . func $
    (-1.5)
  initialise []
  devs <- M.mapM device deviceIDs
  ctxs <- M.mapM (\dev -> CUDA.create dev []) devs
  ptxs <- M.mapM createTargetFromContext ctxs
  let envelope' =
        VU.convert . toUnboxed . computeS $
        analyticalFourierCoefficients2 numR2Freq 1 phiFreq rhoFreq (-s) periodR2
      envelope =
        VS.map
          (\x -> x - VS.sum envelope' / (fromIntegral $ numR2Freq ^ 2))
          envelope'
      envelopeMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        envelope
  plotImageRepaComplex (folderPath </> "EnvelopeFreq.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    envelope
  envelopeR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      numBatchR2
      envelopeMat
  plotImageRepaComplex (folderPath </> "Envelope.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    envelopeR2
      -- envelope2 = VS.map (^ 2) envelope
  let envelope2' =
        VS.map (^ 2) . VU.convert . toUnboxed . computeS $
        analyticalFourierCoefficients2 numR2Freq 1 phiFreq rhoFreq (-s) periodR2
      envelope2 =
        VS.map
          (\x -> x - VS.sum envelope2' / (fromIntegral $ numR2Freq ^ 2))
          envelope2'
      envelope2Mat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        envelope2
  plotImageRepaComplex (folderPath </> "Envelope2Freq.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    envelope2
  envelope2R2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      numBatchR2
      envelope2Mat
  plotImageRepaComplex (folderPath </> "Envelope2.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    envelope2R2
  plan <-
    makePlan
      folderPath
      emptyPlan
      numPointsRecon
      numPointsRecon
      numR2Freq
      (2 * phiFreq + 1)
      (2 * rhoFreq + 1)
  let r2Freqs = L.map fromIntegral . getListFromNumber $ numR2Freq
      -- bias =
      --   VU.convert . toUnboxed . computeUnboxedS $
      --   analyticalFourierCoefficients2 numR2Freq 1 0 0 stdR2 periodR2
      bias' =
        computeUnboxedS $
        analyticalFourierCoefficients2 numR2Freq 1 0 0 stdR2 periodR2
      bias =
        VU.convert .
        toUnboxed .
        computeUnboxedS .
        R.map (\x -> x - sumAllS bias' / (fromIntegral $ numR2Freq ^ 2)) $
        bias' :: VS.Vector (Complex Double)
      biasMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        bias
  plotImageRepaComplex (folderPath </> "BiasFreq.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    bias
  biasR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      numBatchR2
      biasMat
  plotImageRepaComplex (folderPath </> "Bias.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    biasR2
  let forwardPlanID = DFTPlanID DFT1DG [numR2Freq, numR2Freq] [0, 1]
      backwardPlanID = DFTPlanID IDFT1DG [numR2Freq, numR2Freq] [0, 1]
      forwardPlan = getDFTPlan plan forwardPlanID
      backwardPlan = getDFTPlan plan backwardPlanID
      biasVec' =
        VU.convert .
        toUnboxed .
        computeUnboxedS .
        makeFilter2D . fromUnboxed (Z :. numR2Freq :. numR2Freq) . VS.convert $
        bias
      biasVec = VS.map (\x -> x - VS.sum biasVec' / (fromIntegral $ numR2Freq^2)) biasVec'
  biasF <- dftExecuteWithPlan forwardPlanID forwardPlan biasVec
  envelope2F <- dftExecuteWithPlan forwardPlanID forwardPlan envelope2
  biased' <-
    dftExecuteWithPlan backwardPlanID backwardPlan $
    VS.zipWith (*) biasF envelope2F
  let biased =
        VS.map
          (\x -> x - VS.sum biased' / (fromIntegral $ numR2Freq ^ 2))
          biased'
      biasedMat =
        A.use .
        A.fromList (A.Z A.:. (numR2Freq ^ 2) A.:. (1 :: Int)) . VS.toList $
        biased
  plotImageRepaComplex (folderPath </> "BiasedFreq.png") .
    ImageRepa 8 .
    fromUnboxed (Z :. (1 :: Int) :. numR2Freq :. numR2Freq) . VS.convert $
    biased
  biasedR2 <-
    computeFourierSeriesR2StreamAcc
      ptxs
      numR2Freq
      numPointsRecon
      1
      periodR2
      deltaRecon
      numBatchR2
      biasedMat
  plotImageRepaComplex (folderPath </> "Biased.png") .
    ImageRepa 8 .
    computeS . extend (Z :. (1 :: Int) :. All :. All) . sumS . rotate3D $
    biasedR2
