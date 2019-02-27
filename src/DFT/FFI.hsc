{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE ForeignFunctionInterface #-}
module DFT.FFI where

import qualified Foreign.C.Types as C
import Foreign.C.String (CString)
import Foreign.Ptr (Ptr)
import Data.Complex (Complex)

#include <fftw3.h>


type FFTWFlag = C.CUInt

#{enum FFTWFlag,
 , c_measure         = FFTW_MEASURE
 , c_destroy_input   = FFTW_DESTROY_INPUT
 , c_unaligned       = FFTW_UNALIGNED
 , c_conserve_memory = FFTW_CONSERVE_MEMORY
 , c_exhaustive      = FFTW_EXHAUSTIVE
 , c_preserve_input  = FFTW_PRESERVE_INPUT
 , c_patient         = FFTW_PATIENT
 , c_estimate        = FFTW_ESTIMATE
 }


type FFTWSign = C.CInt

#{enum FFTWSign,
 , c_forward = FFTW_FORWARD
 , c_backward = FFTW_BACKWARD
 }


type FFTWKind = C.CInt

#{enum FFTWKind,
 , c_r2hc    = FFTW_R2HC
 , c_hc2r    = FFTW_HC2R
 , c_dht     = FFTW_DHT
 , c_redft00 = FFTW_REDFT00
 , c_redft10 = FFTW_REDFT10
 , c_redft01 = FFTW_REDFT01
 , c_redft11 = FFTW_REDFT11
 , c_rodft00 = FFTW_RODFT00
 , c_rodft10 = FFTW_RODFT10
 , c_rodft01 = FFTW_RODFT01
 , c_rodft11 = FFTW_RODFT11
 }


-- | A plan is an opaque foreign object.
type Plan = Ptr FFTWPlan

type FFTWPlan = ()

-- We use "safe" calls for anything which could take a while so that it won't block
-- other Haskell threads.

-- | Simple plan execution
foreign import ccall safe "fftw3.h fftw_execute" c_execute
    :: Plan -> IO ()
    

-- Execute a plan on different memory than the plan was created for.
-- Alignment /must/ be the same.  If we parallelize a transform of
-- multi-dimensional data by making separate calls within an un-transformed
-- dimension, it is possible that the alignment constraint would not be
-- fulfilled.  However, this only poses a problem for real transforms with odd
-- transform dimension.
foreign import ccall safe "fftw3.h fftw_execute_dft" c_execute_dft
    :: Plan -> Ptr (Complex Double) -> Ptr (Complex Double) -> IO ()
foreign import ccall safe "fftw3.h fftw_execute_dft_r2c" c_execute_dft_r2c
    :: Plan -> Ptr Double -> Ptr (Complex Double) -> IO ()
foreign import ccall safe "fftw3.h fftw_execute_dft_c2r" c_execute_dft_c2r
    :: Plan -> Ptr (Complex Double) -> Ptr Double -> IO ()
foreign import ccall safe "fftw3.h fftw_execute_r2r" c_execute_r2r
    :: Plan -> Ptr Double -> Ptr Double -> IO ()

foreign import ccall safe "fftw3.h fftw_export_wisdom_to_string"
        c_export_wisdom_string :: IO CString

foreign import ccall safe "fftw3.h fftw_import_wisdom_from_string"
        c_import_wisdom_string :: CString -> IO C.CInt

foreign import ccall safe "fftw3.h fftw_import_system_wisdom"
        c_import_wisdom_system :: IO C.CInt

-- | Frees memory allocated by 'fftw_malloc'.  Currently, we only need this to
-- free the wisdom string.
foreign import ccall safe "fftw3.h fftw_free" c_free :: Ptr a -> IO ()


foreign import ccall safe "fftw3.h fftw_destroy_plan" c_destroy_plan ::  Plan -> IO ()

foreign import ccall safe "fftw3.h fftw_cleanup" c_cleanup :: IO ()

foreign import ccall safe "fftw3.h fftw_plan_dft_2d" c_plan_dft_2d :: C.CInt -> C.CInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> C.CInt -> FFTWFlag -> IO Plan

foreign import ccall safe "fftw3.h fftw_plan_dft_1d" c_plan_dft_1d :: C.CInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> C.CInt -> FFTWFlag -> IO Plan

foreign import ccall safe "fftw3.h fftw_plan_dft_r2c_1d" c_plan_dft_r2c_1d :: C.CInt -> Ptr C.CDouble -> Ptr (Complex Double) -> FFTWFlag -> IO Plan

foreign import ccall safe "fftw3.h fftw_plan_many_dft" c_plan_many_dft :: C.CInt -> Ptr C.CInt -> C.CInt -> Ptr (Complex Double) -> Ptr C.CInt -> C.CInt -> C.CInt -> Ptr (Complex Double) -> Ptr C.CInt -> C.CInt -> C.CInt -> C.CInt -> FFTWFlag -> IO Plan

foreign import ccall safe "fftw3.h fftw_plan_dft_r2c_2d" c_plan_dft_r2c_2d :: C.CInt -> C.CInt -> Ptr C.CDouble -> Ptr (Complex Double) -> FFTWFlag -> IO Plan
