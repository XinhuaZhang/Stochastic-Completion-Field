{-# LANGUAGE ViewPatterns #-}
module EigenDFT where

import           Control.Monad                      as M
import           Control.Monad.Parallel             as MP
import qualified Data.Array.Accelerate              as A
import qualified Data.Array.Accelerate.Data.Complex as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Repa                    as R
import           Data.List                          as L
import           Foreign.CUDA.Driver                as CUDA
import           Image.IO
import           Numeric.LinearAlgebra              as NL
import           Pinwheel.Base
import           Pinwheel.FourierSeries2D
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf

main = do
  args@(nStr:numFreqStr:deltaStr:_) <- getArgs
  let n = read nStr :: Int
      numFreq = read numFreqStr :: Int
      delta = read deltaStr :: Double
      len = n ^ 2
      center = div n 2
      centerFreq = div numFreq 2
      constant = (-2) * pi * delta / (fromIntegral n) :: Double
      matArr' =
        fromFunction (Z :. n :. n :. numFreq :. numFreq) $ \(Z :. x' :. y' :. wx' :. wy') ->
          let wx = wx' - centerFreq
              wy = wy' - centerFreq
              x = x' - center
              y = y' - center
          in cis $ constant * fromIntegral (wx * x + wy * y)
      pinwheelVec =
        NL.fromList . R.toList . fromFunction (Z :. n :. n) $ \(Z :. x' :. y') ->
          let x = x' - center
              y = y' - center
          in fourierMellin (0.5 :: Double) 5 5 (fromIntegral x, fromIntegral y)
      folderPath = "output/test/EigenDFT"
  when (numFreq > n) (error "error: numFreq > n")
  removePathForcibly folderPath
  createDirectoryIfMissing True folderPath
  matArr <- computeUnboxedP matArr'
  let mat = (len >< (numFreq ^ 2)) . R.toList $ matArr
      coefficients = pinwheelVec <# mat
  --     xs =
  --       if numFreq == n
  --         then let (eigVal, eigVec) = eig mat
  --              in snd . L.unzip . L.reverse . L.sortOn fst $
  --                 L.zip
  --                   (L.map magnitude . NL.toList $ eigVal)
  --                   (L.map NL.toList . NL.toColumns $ eigVec)
  --         else let (u, s, v) = thinSVD mat
  --              in snd . L.unzip . L.reverse . L.sortOn fst $
  --                 L.zip (NL.toList $ s) (L.map NL.toList . NL.toColumns $ u)
  -- M.mapM_
  --   (\(i, ys) ->
  --      plotImageRepaComplex (folderPath </> printf "%d.png" i) .
  --      ImageRepa 8 . fromListUnboxed (Z :. (1 :: Int) :. n :. n) $
  --      ys) $
  --   L.zip [0 :: Int .. (numFreq ^ 2 - 1)] xs
  plotImageRepaComplex (folderPath </> "Pinwheel.png" ) . ImageRepa 8 . fromListUnboxed (Z :. (1 :: Int) :. n :. n) . NL.toList $ pinwheelVec
  plotImageRepaComplex (folderPath </> "coefficients.png" ) . ImageRepa 8 . fromListUnboxed (Z :. (1 :: Int) :. numFreq :. numFreq) . NL.toList $ coefficients
  -- let centerAcc = A.constant center
  --     centerFreqAcc = A.constant centerFreq
  --     constantAcc = A.constant constant
  --     matArrAcc' =
  --       A.generate (A.constant (A.Z A.:. numFreq A.:. numFreq A.:. n A.:. n)) $ \(A.unlift -> (A.Z A.:. wx' A.:. wy' A.:. x' A.:. y')) ->
  --         let wx = wx' - centerFreqAcc
  --             wy = wy' - centerFreqAcc
  --             x = x' - centerAcc
  --             y = y' - centerAcc
  --         in A.cis $ constantAcc * A.fromIntegral (wx * x + wy * y)
  -- initialise []
  -- dev <- device 7
  -- ctx <- CUDA.create dev []
  -- ptx <- createTargetFromContext ctx
  -- let matArrAcc = runWith ptx matArrAcc'
  --     mat = ((numFreq^2) >< len) . A.toList $ matArrAcc
