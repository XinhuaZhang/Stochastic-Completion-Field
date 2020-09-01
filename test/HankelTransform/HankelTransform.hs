module HankelTransform where

import           Control.Parallel.Strategies
import           Data.Array.Repa                        as R
import           Data.Complex                           as C
import           Data.List                              as L
import           Data.Vector.Storable                   as VS
import           Data.Vector.Unboxed                    as VU
import           Graphics.Rendering.Chart.Backend.Cairo
import           Graphics.Rendering.Chart.Easy
import           Numeric.GSL.Special.Bessel
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.BLAS
import           Utils.List
import           Utils.Parallel

main = do
  args@(nStr:deltaXStr:deltaRStr:radiusXStr:radiusRStr:_) <- getArgs
  let n = read nStr :: Int
      deltaX = read deltaXStr :: Double
      deltaR = read deltaRStr :: Double
      radiusX = read radiusXStr :: Double
      radiusR = read radiusRStr :: Double
      folderPath = "output/test/HankelTransform"
  createDirectoryIfMissing True folderPath
  let numR = round $ radiusR / deltaR :: Int
      numX = round $ radiusX / deltaX :: Int
      func r =
        if r >= 1
          then 1 / r
          else 1 / r
      fr =
        parMap
          rdeepseq
          (\r -> func r * r)
          [fromIntegral r * deltaR | r <- [1 .. numR]]
      besselArrayR =
        fromFunction (Z :. numX :. numR) $ \(Z :. x :. r) ->
          deltaR *
          bessel_Jn
            (fromIntegral n)
            ((fromIntegral x + 1) * deltaX * (fromIntegral r + 1) * deltaR)
  besselArray <- computeUnboxedP besselArrayR
  let besselVec = VU.convert . toUnboxed $ besselArray
      frVec = VS.fromList fr
  xs <- VS.toList <$> gemmBLAS numX 1 numR besselVec frVec
  toFile def (folderPath </> "HankelTransform.png") $ do
    layout_title .=
      printf "N = %d RadiusR = %.2f RadiusX = %.2f" n radiusR radiusX
    plot (line "" [L.zip [fromIntegral i * deltaX | i <- [1 .. numX]] xs])
