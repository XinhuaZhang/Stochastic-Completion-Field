{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE Strict            #-}
{-# LANGUAGE StrictData        #-}
module Utils.BLAS where

import           Data.Complex
import           Data.Vector.Storable            as VS
import           Foreign.CUDA.BLAS               as BLAS
import           Foreign.CUDA.Driver             as CUDA
import           Foreign.Marshal
import           Numerical.HBLAS.BLAS.FFI
import           Numerical.HBLAS.BLAS.FFI.Level1
import           Numerical.HBLAS.BLAS.FFI.Level3

type BLASMMT e
   = Int -> Int -> Int -> VS.Vector e -> VS.Vector e -> IO (VS.Vector e)
   
type BLASDOTU e = VS.Vector e -> VS.Vector e -> IO e

class BLAS a where
  gemmBLAS :: BLASMMT a
  dotuBLAS :: BLASDOTU a

 
instance BLAS Double where
  {-# INLINE gemmBLAS #-}
  gemmBLAS m' n' k' a b = do
    let output = VS.replicate (m' * n') 0
        !m = fromIntegral m'
        !n = fromIntegral n'
        !k = fromIntegral k'
    unsafeWith a $ \aPtr ->
      unsafeWith b $ \bPtr ->
        unsafeWith output $ \cPtr ->
          cblas_dgemm_safe
            (encodeOrder BLASRowMajor)
            (encodeTranspose BlasNoTranspose)
            (encodeTranspose BlasNoTranspose)
            m
            n
            k
            1
            aPtr
            k
            bPtr
            n
            0
            cPtr
            n
    return output
  {-# INLINE dotuBLAS #-}
  dotuBLAS vec1 vec2 =
    unsafeWith vec1 $ \ptr1 ->
      unsafeWith vec2 $ \ptr2 ->
        cblas_ddot_safe (fromIntegral . VS.length $ vec1) ptr1 1 ptr2 1

instance BLAS (Complex Double) where
  {-# INLINE gemmBLAS #-}
  gemmBLAS m' n' k' a b = do
    let output = VS.replicate (m' * n') 0
        !m = fromIntegral m'
        !n = fromIntegral n'
        !k = fromIntegral k'
    with 1 $ \alphaPtr ->
      with 0 $ \betaPtr ->
        unsafeWith a $ \aPtr ->
          unsafeWith b $ \bPtr ->
            unsafeWith output $ \cPtr ->
              cblas_zgemm_safe
                (encodeOrder BLASRowMajor)
                (encodeTranspose BlasNoTranspose)
                (encodeTranspose BlasNoTranspose)
                m
                n
                k
                alphaPtr
                aPtr
                k
                bPtr
                n
                betaPtr
                cPtr
                n
    return output
  {-# INLINE dotuBLAS #-}
  dotuBLAS vec1 vec2 = do
    let output = VS.singleton 0
    unsafeWith vec1 $ \ptr1 ->
      unsafeWith vec2 $ \ptr2 ->
        unsafeWith output $ \outPtr ->
          cblas_zdotu_safe
            (fromIntegral . VS.length $ vec1)
            ptr1
            1
            ptr2
            1
            outPtr
    return . VS.head $ output



type CUBLASMMT a b e
   = Handle -> Int -> Int -> Int -> a e -> b e -> IO (VS.Vector e)

class CUBLAS a where
  gemmCuBLAS :: CUBLASMMT DevicePtr DevicePtr a

instance CUBLAS Float where
  {-# INLINE gemmCuBLAS #-}
  gemmCuBLAS handle !m !n !k matA matB = do
    let !sizeC = m * n
        output = VS.replicate sizeC 0
    with 1 $ \alpha ->
      with 0 $ \beta ->
        CUDA.allocaArray sizeC $ \matC ->
          unsafeWith output $ \outputPtr -> do
            sgemm handle N N n m k alpha matB n matA k beta matC n
            CUDA.peekArray sizeC matC outputPtr
    return output

instance CUBLAS Double where
  {-# INLINE gemmCuBLAS #-}
  gemmCuBLAS handle !m !n !k matA matB = do
    let !sizeC = m * n
        output = VS.replicate sizeC 0
    with 1 $ \alpha ->
      with 0 $ \beta ->
        CUDA.allocaArray sizeC $ \matC ->
          unsafeWith output $ \outputPtr -> do
            dgemm handle N N n m k alpha matB n matA k beta matC n
            CUDA.peekArray sizeC matC outputPtr
    return output

instance CUBLAS (Complex Float) where
  {-# INLINE gemmCuBLAS #-}
  gemmCuBLAS handle !m !n !k matA matB = do
    let !sizeC = m * n
        output = VS.replicate sizeC 0
    with 1 $ \alpha ->
      with 0 $ \beta ->
        CUDA.allocaArray sizeC $ \matC ->
          unsafeWith output $ \outputPtr -> do
            cgemm handle N N n m k alpha matB n matA k beta matC n
            CUDA.peekArray sizeC matC outputPtr
    return output

instance CUBLAS (Complex Double) where
  {-# INLINE gemmCuBLAS #-}
  gemmCuBLAS handle !m !n !k matA matB = do
    let !sizeC = m * n
        output = VS.replicate sizeC 0
    with 1 $ \alpha ->
      with 0 $ \beta ->
        CUDA.allocaArray sizeC $ \matC ->
          unsafeWith output $ \outputPtr -> do
            zgemm handle N N n m k alpha matB n matA k beta matC n
            CUDA.peekArray sizeC matC outputPtr
    return output

{-# INLINE unsafeWithGPU #-}
unsafeWithGPU :: (Storable e) => VS.Vector e -> (DevicePtr e -> IO b) -> IO b
unsafeWithGPU vec f = do
  let !len = VS.length vec
  CUDA.allocaArray len $ \devPtr -> do
    unsafeWith vec $ \ptr -> CUDA.pokeArray len ptr devPtr
    f devPtr

{-# INLINE gemmCuBLAS10 #-}
gemmCuBLAS10 :: (CUBLAS e, Storable e) => CUBLASMMT VS.Vector DevicePtr e
gemmCuBLAS10 handle m n k a matB =
  unsafeWithGPU a $ \matA -> gemmCuBLAS handle m n k matA matB

{-# INLINE gemmCuBLAS01 #-}
gemmCuBLAS01 :: (CUBLAS e, Storable e) => CUBLASMMT DevicePtr VS.Vector e
gemmCuBLAS01 handle m n k matA b =
  unsafeWithGPU b $ \matB -> gemmCuBLAS handle m n k matA matB

{-# INLINE gemmCuBLAS11 #-}
gemmCuBLAS11 :: (CUBLAS e, Storable e) => CUBLASMMT VS.Vector VS.Vector e
gemmCuBLAS11 handle m n k a b =
  unsafeWithGPU a $ \matA ->
    unsafeWithGPU b $ \matB -> gemmCuBLAS handle m n k matA matB
