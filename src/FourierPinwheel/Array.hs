{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE StrictData       #-}
module FourierPinwheel.Array where

import qualified Data.Array.Accelerate as A
import Data.Array.Repa       as R
import           Data.Complex
import           Data.List             as L
import           Data.Vector.Generic   as VG
import           Data.Vector.Storable  as VS
import           DFT.Plan
import           Filter.Utils
import           Image.IO
import           Utils.Array
import           Utils.Parallel
import           Graphics.Rendering.Chart.Backend.Cairo
import           Graphics.Rendering.Chart.Easy
import           Utils.BLAS
import Text.Printf

data FPArray vector = FPArray
  { getFPArrayNumXFreq     :: Int
  , getFPArrayNumYFreq     :: Int
  , getFPArrayNumRFreq     :: Int
  , getFPArrayNumThetaFreq :: Int
  , getFPArrayNumRhoFreq   :: Int
  , getFPArrayNumPhiFreq   :: Int
  , getFPArray             :: [vector]
  }

{-# INLINE toMatrixAcc #-}
toMatrixAcc ::
     (VG.Vector vector e, A.Elt e) => FPArray (vector e) -> A.Acc (A.Array A.DIM2 e)
toMatrixAcc (FPArray numXFreq numYFreq numRFreq numThetaFreq _ _ arr) =
  A.use .
  A.fromList (A.Z A.:. (numRFreq * numThetaFreq) A.:. (numXFreq * numYFreq)) .
  VG.toList . VG.concat $
  arr

{-# INLINE parMapFPArray #-}
parMapFPArray ::
     (NFData vector) => (vector -> vector) -> FPArray vector -> FPArray vector
parMapFPArray f (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs) =
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
  parMap rdeepseq f $
  vecs

{-# INLINE parZipWithFPArray #-}
parZipWithFPArray ::
     (NFData vector)
  => (vector -> vector -> vector)
  -> FPArray vector
  -> FPArray vector
  -> FPArray vector
parZipWithFPArray f (FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq vecs1) =
  FPArray numXFreq numYFreq numRFreq numThetaFreq numRhoFreq numPhiFreq .
  parZipWith rdeepseq f vecs1 . getFPArray 


plotFPArray ::
     DFTPlan
  -> FilePath
  -> FPArray (VS.Vector (Complex Double))
  -> IO (R.Array R.U R.DIM4 (Complex Double)) 
plotFPArray plan filePath arr = do
  let planIDBackward =
        DFTPlanID
          IDFT1DG
          [ getFPArrayNumThetaFreq arr
          , getFPArrayNumXFreq arr
          , getFPArrayNumYFreq arr
          ]
          [1, 2]
  vecsR2 <- dftExecuteBatchP plan planIDBackward . getFPArray $ arr
  arrR2 <-
    R.computeUnboxedP .
    makeFilter2DInverse .
    R.fromUnboxed
      (R.Z R.:. getFPArrayNumRFreq arr R.:. getFPArrayNumThetaFreq arr R.:.
       getFPArrayNumXFreq arr R.:.
       getFPArrayNumYFreq arr) .
    VS.convert . VS.concat $
    vecsR2
  arr1 <-
    R.sumP . R.sumS . R.map (\x -> magnitude x Prelude.^ 2) . rotate4D2 $ arrR2
  plotImageRepa filePath .
    ImageRepa 8 .
    R.fromUnboxed
      (R.Z R.:. (1 :: Int) R.:. getFPArrayNumXFreq arr R.:.
       getFPArrayNumYFreq arr) .
    -- VG.map (\x -> x^2) . 
    R.toUnboxed $ arr1
  return arrR2


plotRThetaDist ::
     FilePath
  -> FilePath
  -> Int
  -> Int
  -> Int
  -> Double
  -> (Int, Int)
  -> R.Array R.U R.DIM4 (Complex Double)
  -> IO ()
plotRThetaDist filePathTheta filePathR numPoints numTheta numR periodEnv (x, y) arr = do
  let (Z :. numRFreq :. numThetaFreq :. _ :. _) = extent arr
      center = div numPoints 2
      centerR = div numR 2
      centerTheta = div numTheta 2
      centerRFreq = div numRFreq 2
      centerThetaFreq = div numThetaFreq 2
      vecRTheta =
        VG.convert . R.toUnboxed . R.computeUnboxedS . R.slice arr $
        (R.Z R.:. R.All R.:. R.All R.:. (x + center) R.:. (y + center))
      deltaTheta = 2 * pi / Prelude.fromIntegral numTheta
      deltaLogR = log periodEnv / Prelude.fromIntegral numR
      matrix =
        VG.convert .
        R.toUnboxed .
        R.computeUnboxedS .
        R.fromFunction (R.Z R.:. numR R.:. numTheta :. numRFreq :. numThetaFreq) $ \(Z :. r :. t :. rf :. tf) ->
          cis
            (fromIntegral ((r - centerR) * (rf - centerRFreq)) * 2 * pi /
             log periodEnv *
             deltaLogR +
             fromIntegral ((t - centerTheta) * (tf - centerThetaFreq)) *
             deltaTheta)
      m = numR * numTheta
      n = 1
      k = numRFreq * numThetaFreq
  vec <- VS.map magnitude <$> gemmBLAS m n k matrix vecRTheta
  let outputArr = fromUnboxed (Z :. numR :. numTheta) . VG.convert $ vec
      outputR = R.toList . sumS $ outputArr
      outputTheta =
        R.toList .
        sumS .
        R.backpermute (Z :. numTheta :. numR) (\(Z :. i :. j) -> (Z :. j :. i)) $
        outputArr
  toFile def filePathTheta $ do
    layout_title .= printf "Theta (%d , %d)" x y
    plot
      (line
         ""
         [ L.zip
             [fromIntegral (i - centerTheta) * deltaTheta | i <- [0 .. numTheta - 1]]
             outputTheta
         ])
  toFile def filePathR $ do
    layout_title .= printf "LogR (%d , %d)" x y
    plot
      (line
         ""
         [ L.zip
             [fromIntegral (i - centerR) * deltaLogR | i <- [0 .. numR - 1]]
             outputR
         ])
