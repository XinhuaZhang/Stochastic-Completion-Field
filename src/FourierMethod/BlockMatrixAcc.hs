{-# LANGUAGE Strict #-}
module FourierMethod.BlockMatrixAcc where

import           Data.Array.Accelerate                        as A
import           Data.Array.Accelerate.LLVM.PTX
import           Data.Array.Accelerate.Numeric.LinearAlgebra  as A
import           Data.List                                    as L
import           Data.Vector.Generic                          as VG
import           Data.Vector.Unboxed                          as VU
import           Text.Printf
import           Utils.Parallel

{-# INLINE getRows #-}
getRows :: Matrix e -> Int
getRows arr =
  let (Z :. rows :. _) = arrayShape arr
  in rows

{-# INLINE getCols #-}
getCols :: Matrix e -> Int
getCols arr =
  let (Z :. _ :. cols) = arrayShape arr
  in cols

{-# INLINE blockMatrixMultiply #-}
blockMatrixMultiply ::
     (Numeric e, Unbox e)
  => [PTX]
  -> [Acc (Matrix e)]
  -> Acc (Matrix e)
  -> VU.Vector e
blockMatrixMultiply ptxs matAs matB =
  if (L.length ptxs Prelude.== L.length matAs)
    then VU.concat .
         parZipWith
           rdeepseq
           (\ptx matA -> VU.fromList . A.toList $ runWith ptx (matA A.<> matB))
           ptxs $
         matAs
    else error $
         printf
           "Error in blockMatrixMultiply: matrix A is not divided according to the number of GPUs.\n%d GPUs vs %d blocks"
           (L.length ptxs)
           (L.length matAs)

{-# INLINE fromDFTArray #-}
fromDFTArray :: (VG.Vector vector e, Elt e) => [vector e] -> Matrix e
fromDFTArray xs =
  let rows = L.length xs
      cols = VG.length . L.head $ xs
  in A.fromList (Z :. rows :. cols) . VG.toList . VG.concat $ xs
