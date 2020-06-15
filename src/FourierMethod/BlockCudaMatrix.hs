{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
module FourierMethod.BlockCudaMatrix where

import           Control.Concurrent.Async
import           Control.DeepSeq
import           Control.Monad               as M
import           Control.Monad.Parallel               as MP
import           Data.Array.Repa             as R
import           Data.Complex
import           Data.List                   as L
import           Data.Vector.Storable        as VS
import           Data.Vector.Unboxed         as VU
import           Foreign.CUDA.BLAS           as BLAS
import           Foreign.CUDA.Driver         as CUDA
import           Foreign.CUDA.Runtime.Device as CudaRT
import           Text.Printf
import           Utils.BLAS
import           Utils.Parallel
import Control.Concurrent.MVar

data CuVec a
  = CuVecHost (VS.Vector a)
  | CuVecDevice (DevicePtr a)

data CuMat a =
  CuMat Int
        Int
        (CuVec a)

instance (Storable a, Show a) => Show (CuVec a) where
  show (CuVecHost vec) = printf "CuVecHost length %d %s" (VS.length vec) (show . VS.head $ vec)
  show (CuVecDevice _) = "CuVecDevice"

instance NFData (CuVec a) where
  rnf (CuVecHost vec) = seq vec ()
  rnf (CuVecDevice ptr) = seq ptr ()

instance (Storable a, Show a) => Show (CuMat a) where
  show (CuMat rows cols vec) =
    printf "CuMat rows = %d cols = %d %s" rows cols (show vec)

instance NFData (CuMat a) where
  rnf (CuMat rows cols vec) = rows `seq` cols `seq` vec `seq` ()

{-# INLINE getHostCuVec #-}
getHostCuVec :: (Storable a) => CuVec a -> VS.Vector a
getHostCuVec (CuVecHost vec) = vec
getHostCuVec (CuVecDevice _) = error "Error in getHostCuVec: no host vector."

{-# INLINE getHostCuMat #-}
getHostCuMat :: (Storable a) => CuMat a -> VS.Vector a
getHostCuMat (CuMat _ _ vec) = getHostCuVec vec

{-# INLINE getRowsCuMat #-}
getRowsCuMat :: CuMat a -> Int
getRowsCuMat (CuMat rows _ _) = rows

{-# INLINE getColsCuMat #-}
getColsCuMat :: CuMat a -> Int
getColsCuMat (CuMat _ cols _) = cols

{-# INLINE concatCuMat #-}
concatCuMat :: (Storable a) => [CuMat a] -> CuMat a
concatCuMat xs =
  let rows = L.foldl' (\s mat -> s + getRowsCuMat mat) 0 xs
      cols = getColsCuMat . L.head $ xs
  in if L.any (/= cols) . L.map getColsCuMat $ xs
       then error $
            printf
              "concatCuMat: concat matries with different numbers of columns.\n%s\n"
              (show . L.map getColsCuMat $ xs)
       else CuMat rows cols . CuVecHost . VS.concat . L.map getHostCuMat $ xs

{-# INLINE unsafeWithCuMat #-}
unsafeWithCuMat :: (Storable e) => CuMat e -> (DevicePtr e -> IO a) -> IO a
unsafeWithCuMat (CuMat _ _ (CuVecHost vec)) f = unsafeWithGPU vec f
unsafeWithCuMat (CuMat _ _ (CuVecDevice vec)) f = f vec

{-# INLINE freeCuMat #-}
freeCuMat :: CuMat e -> IO ()
freeCuMat (CuMat _ _ (CuVecHost _)) = return ()
freeCuMat (CuMat _ _ (CuVecDevice ptr)) = CUDA.free ptr

{-# INLINE subMatMul #-}
subMatMul ::
     (Storable e, Floating e, CUBLAS e, Show e)
  => Handle
  -> Int
  -> [CuMat e]
  -> CuMat e
  -> IO (CuMat e)
subMatMul handle rows matA@((CuMat rowsA colsA _):_) b@(CuMat _ colsB _) =
  fmap (CuMat rows colsB . CuVecHost . VS.concat) . unsafeWithCuMat b $ \bPtr ->
    M.mapM
      (\harmonic ->
         unsafeWithCuMat harmonic $ \harmonicPtr ->
           gemmCuBLAS handle rowsA colsB colsA harmonicPtr bPtr)
      matA

-- Given row-major Mat1: rows1xcols1 and Mat2: rows2xcols2, this function outputs Mat1 X Mat2
-- {-# INLINE matMul #-}
matMul ::
     (Storable e, Floating e, CUBLAS e, Unbox e, Show e)
  => MVar ()
  -> Bool
  -> Int
  -> Int
  -> [CuMat e]
  -> [CuMat e]
  -> IO (CuMat e)
matMul lock doesTranspose deviceID rows matA' matB' = do
  a <- takeMVar lock
  CudaRT.set deviceID
  handle <- BLAS.create
  putMVar lock a
  matA <-
    if L.length matA' == 1
      then do
        let (CuMat rowsA colsA (CuVecHost vec)) = L.head matA'
            len = rowsA * colsA
        devPtr <- CUDA.mallocArray len
        unsafeWith vec $ \ptr -> CUDA.pokeArray len ptr devPtr
        return [CuMat rowsA colsA (CuVecDevice devPtr)]
      else return matA'
  matB <-
    if L.length matB' == 1
      then do
        let (CuMat rowsB colsB (CuVecHost vec)) = L.head matB'
            len = rowsB * colsB
        devPtr <- CUDA.mallocArray len
        unsafeWith vec $ \ptr -> CUDA.pokeArray len ptr devPtr
        return [CuMat rowsB colsB (CuVecDevice devPtr)]
      else return matB'
  coefs <- M.mapM (subMatMul handle rows matA) matB
  let transposedCoefs = concatCuMat . L.map transposeCuMat $!! coefs
      output =
        if doesTranspose
          then transposedCoefs
          else transposeCuMat transposedCoefs
  when (L.length matA' == 1) (freeCuMat . L.head $ matA)
  when (L.length matB' == 1) (freeCuMat . L.head $ matB)
  BLAS.destroy handle
  return output

-- matA is divided into [rowsA1 x colsA .. rowsAN x colsA]
-- matB is divided into [colsA x colsB1 .. colsA x colsBM]
-- each colsA x colsBm is further divided into [colsA x colsBmk]
-- M is the number of GPUs
{-# INLINE blockMatrixMultiply1 #-}
blockMatrixMultiply1 ::
     (Storable e, Floating e, CUBLAS e, Unbox e, Show e)
  => Bool
  -> [Int]
  -> [CuMat e]
  -> [[CuMat e]]
  -> IO (CuMat e)
blockMatrixMultiply1 doesTranspose deviceIDs matAs matBss = do
  unless
    (L.length deviceIDs == L.length matBss)
    (error $
     printf
       "Error in blockMatrixMultiply1: matrix B is not divided according to the number of GPUs.\n%d GPUs vs %d blocks"
       (L.length deviceIDs)
       (L.length matBss))
  let rows = L.sum . L.map getRowsCuMat $ matAs
  lock <- newMVar ()
  output <-
    fmap concatCuMat .
    mapConcurrently (\(deviceID, matBs) -> matMul lock True deviceID rows matAs matBs) .
    L.zip deviceIDs $
    matBss
  return $!!
    if doesTranspose
      then output
      else transposeCuMat output


{-# INLINE blockMatrixMultiply2 #-}
blockMatrixMultiply2 ::
     (Storable e, Floating e, CUBLAS e, Unbox e, Show e)
  => Bool
  -> [Int]
  -> [[CuMat e]]
  -> [CuMat e]
  -> IO (CuMat e)
blockMatrixMultiply2 doesTranspose deviceIDs matAss matBs = do
  unless
    (L.length deviceIDs == L.length matAss)
    (error $
     printf
       "Error in blockMatrixMultiply2: matrix A is not divided according to the number of GPUs.\n%d GPUs vs %d blocks"
       (L.length deviceIDs)
       (L.length matAss))
  let rowss = L.map (L.sum . L.map getRowsCuMat) $ matAss
  lock <- newMVar ()
  output <-
    fmap concatCuMat .
    mapConcurrently
      (\(deviceID, rows, matAs) -> matMul lock False deviceID rows matAs matBs) .
    L.zip3 deviceIDs rowss $
    matAss
  return $!!
    if doesTranspose
      then transposeCuMat output
      else output
      

{-# INLINE blockMatrixMultiply3 #-}
blockMatrixMultiply3 ::
     (Storable e, Floating e, CUBLAS e, Unbox e, Show e)
  => Bool
  -> [Int]
  -> [CuMat e]
  -> [CuMat e]
  -> IO (CuMat e)
blockMatrixMultiply3 doesTranspose deviceIDs matAs matBs = do
  unless
    (L.length deviceIDs == L.length matAs)
    (error $
     printf
       "Error in blockMatrixMultiply3: matrix A is not divided according to the number of GPUs.\n%d GPUs vs %d blocks"
       (L.length deviceIDs)
       (L.length matAs))
  let rowss = L.map getRowsCuMat $ matAs
  lock <- newMVar ()
  output <-
    fmap concatCuMat .
    mapConcurrently (\(deviceID, rows, matA) -> matMul lock False deviceID rows [matA] matBs) .
    L.zip3 deviceIDs rowss $
    matAs
  return $!!
    if doesTranspose
      then transposeCuMat output 
      else output


-- Utilities
{-# INLINE transposeCuMat #-}
transposeCuMat :: (Storable e, Unbox e) => CuMat e -> CuMat e
transposeCuMat (CuMat rows cols (CuVecHost vec)) =
  CuMat cols rows .
  CuVecHost .
  VS.convert .
  toUnboxed .
  computeS .
  R.backpermute (Z :. cols :. rows) (\(Z :. r :. c) -> (Z :. c :. r)) .
  fromUnboxed (Z :. rows :. cols) . VS.convert $
  vec
transposeCuMat _ = error "transposeCuMat: Cannot transpose matrices on a device."

{-# INLINE createCuMat #-}
createCuMat :: (Storable e) => [VS.Vector e] -> CuMat e
createCuMat xs =
  CuMat (L.length xs) (VS.length . L.head $ xs) . CuVecHost . VS.concat $ xs
