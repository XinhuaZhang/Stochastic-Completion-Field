{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE StrictData       #-}
module Pinwheel.Transform where

import           Control.DeepSeq
import           Data.Array.IArray          as IA
import           Data.Complex
import           Data.List                  as L
import           Data.Vector.Generic        as VG
import           Data.Vector.Unboxed        as VU
import           FokkerPlanck.FourierSeries (getHarmonics)
import           FokkerPlanck.Histogram
import           Pinwheel.List
import           Sparse.Vector              as SV
import           Utils.List
import           Utils.Parallel             hiding (dot)

data PinwheelTransformData a b =
  PinwheelTransformData a
                        a
                        b

instance (NFData a, NFData b) => NFData (PinwheelTransformData a b) where
  rnf (PinwheelTransformData x y z) = x `seq` y `seq` z `seq` ()

{-# INLINE projectVec #-}
projectVec ::
     ( VG.Vector vector e
     , VG.Vector vector (Complex e)
     , VG.Vector vector Int
     , Num e
     , RealFloat e
     )
  => Int
  -> Int
  -> e
  -> vector (Complex e)
  -> PinwheelTransformData e (SparseVector vector (Complex e))
  -> Complex e
projectVec angularFreq radialFreq sigma pinwheel (PinwheelTransformData logR theta sparseVec) =
  ((exp ((sigma - 1) * logR)) :+ 0) *
  (cis $
   (-1) * (theta * fromIntegral angularFreq + logR * fromIntegral radialFreq)) *
  (sparseVec `dot` pinwheel)

{-# INLINE pinwheelTransform #-}
pinwheelTransform ::
     ( VG.Vector vector e
     , VG.Vector vector (Complex e)
     , VG.Vector vector Int
     , NFData e
     , RealFloat e
     , Unbox e
     )
  => Int
  -> Int
  -> Int
  -> Int
  -> Int
  -> e
  -> IA.Array (Int, Int) (vector (Complex e))
  -> [PinwheelTransformData e (SparseVector vector (Complex e))]
  -> Histogram (Complex e)
pinwheelTransform numThread maxPhiFreq maxRhoFreq maxThetaFreq maxRFreq sigma pinwheelArr xs =
  let coef =
        parMapChunk
          (ParallelParams numThread undefined)
          rdeepseq
          (\(rFreq, thetaFreq, rhoFreq, phiFreq) ->
             let pinwheel =
                   getHarmonics
                     pinwheelArr
                     (fromIntegral phiFreq)
                     (fromIntegral rhoFreq)
                     (fromIntegral thetaFreq)
                     (fromIntegral rFreq)
             in L.foldl'
                  (\s x -> s + projectVec thetaFreq rFreq sigma pinwheel x)
                  0
                  xs) $
        [ (rFreq, thetaFreq, rhoFreq, phiFreq)
        | rFreq <- [-maxRFreq .. maxRFreq]
        , thetaFreq <- [-maxThetaFreq .. maxThetaFreq]
        , rhoFreq <- [-maxRhoFreq .. maxRhoFreq]
        , phiFreq <- [-maxPhiFreq .. maxPhiFreq]
        ]
  in Histogram
       [ 2 * maxPhiFreq + 1
       , 2 * maxRhoFreq + 1
       , 2 * maxThetaFreq + 1
       , 2 * maxRFreq + 1
       ]
       1 .
     VU.fromList $
     coef
