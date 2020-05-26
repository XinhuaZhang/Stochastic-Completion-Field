{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
module STC.Convolution where

import           Control.Arrow
import           Control.Monad              as M
import           Control.Monad.Parallel     as MP (bindM2, mapM)
import           Data.Array.IArray          as IA
import           Data.Array.Repa            as R
import           Data.Complex
import           Data.Ix
import           Data.List                  as L
import           Data.Vector.Storable       as VS
import           Data.Vector.Unboxed        as VU
import           Debug.Trace
import           DFT.Plan
import           Filter.Utils
import           FokkerPlanck.FourierSeries
import           Image.IO
import           STC.DFTArray
import           Types
import           Utils.Parallel
import FokkerPlanck.GPUKernel
import qualified Data.Array.Accelerate              as A
import qualified Data.Array.Accelerate.Data.Complex as A
import           Data.Array.Accelerate.LLVM.PTX
import Debug.Trace

data Field
  = Source
  | Sink
  deriving (Read, Show)

{-# INLINE dftHarmonicsArray #-}
dftHarmonicsArray ::
     DFTPlan
  -> Int
  -> Double
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double
  -> IO (IA.Array (Int, Int) (VS.Vector (Complex Double)))
dftHarmonicsArray !plan !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff = do
  let !arr =
        computeHarmonicsArray
          numRows
          deltaRow
          numCols
          deltaCol
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          halfLogPeriod
          cutoff
  fmap (listArray (bounds arr)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [numCols, numRows] [0, 1]) .
    L.map
      (VS.convert .
       toUnboxed .
       computeUnboxedS . makeFilter2D . fromUnboxed (Z :. numCols :. numRows)) .
    IA.elems $
    arr
  
{-# INLINE dftHarmonicsArrayG #-}
dftHarmonicsArrayG ::
     DFTPlan
  -> Int
  -> Double
  -> Int
  -> Double
  -> [Double]
  -> [Double]
  -> [Double]
  -> [Double]
  -> Double
  -> Double -> VS.Vector (Complex Double)
  -> IO (IA.Array (Int, Int) (VS.Vector (Complex Double)))
dftHarmonicsArrayG !plan !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff !gaussian = do
  let !arr =
        computeHarmonicsArray
          numRows
          deltaRow
          numCols
          deltaCol
          phiFreqs
          rhoFreqs
          thetaFreqs
          rFreqs
          halfLogPeriod
          cutoff
  xs <-
    fmap (L.map (VS.zipWith (*) gaussian)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [numCols, numRows] [0, 1]) .
    IA.elems $
    arr
  ys <- dftExecuteBatchP plan (DFTPlanID IDFT1DG [numCols, numRows] [0, 1]) xs
  fmap (listArray (bounds arr)) .
    dftExecuteBatchP plan (DFTPlanID DFT1DG [numCols, numRows] [0, 1]) .
    L.map
      (\vec ->
         let arr = fromUnboxed (Z :. numCols :. numRows) . VS.convert $ vec
         in VS.convert .
            toUnboxed . computeUnboxedS . makeFilter2D . R.traverse arr id $ \f idx@(Z :. i :. j) ->
              if (i == div numCols 2 && j == div numRows 2)
                 then 0
                 else f idx) $
    ys


{-# INLINE convolve #-}
convolve ::
     Field
  -> DFTPlan 
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> DFTArray
  -> IO DFTArray
convolve !field !plan !coefficients !harmonicsArray !arr@(DFTArray rows cols thetaFreqs rFreqs vecs) = do
  dftVecs <- dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) vecs
  let !initVec = VS.replicate (VS.length . L.head $ vecs) 0
      idx = (,) <$> (L.zip [0 ..] rFreqs) <*> (L.zip [0 ..] thetaFreqs)
  fmap (DFTArray rows cols thetaFreqs rFreqs) .
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) .
    parMap
      rdeepseq
      (\((!r, !rFreq), (!theta, !thetaFreq)) ->
         VS.map
           (\x ->
              case field of
                Source -> x
                Sink -> x * cis (thetaFreq * pi)) .
         L.foldl'
           (\(!vec) (((!rho, !rhoFreq), (!phi, !phiFreq)), inputVec) ->
              VS.zipWith
                (+)
                vec
                (VS.map (* (coefficients R.! (Z :. r :. theta :. rho :. phi))) .
                 VS.zipWith
                   (*)
                   (getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq rFreq) $
                 inputVec))
           initVec .
         L.zip idx $
         dftVecs) $
    idx

{-# INLINE convolve' #-}
convolve' ::
     Field
  -> DFTPlan
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> DFTArray
  -> IO DFTArray
convolve' !field !plan !coefficients !harmonicsArray !arr@(DFTArray rows cols thetaFreqs rFreqs vecs) = do
  dftVecs <- dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) vecs
  let !initVec = VS.replicate (VS.length . L.head $ vecs) 0
      idxTheta = L.zip [0 ..] thetaFreqs
  fmap (DFTArray rows cols thetaFreqs rFreqs) .
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) .
    parMap
      rdeepseq
      (\(!theta, !thetaFreq) ->
         VS.map
           (\x ->
              case field of
                Source -> x
                Sink -> x * cis (pi * thetaFreq)) .
         L.foldl'
           (\vec1 (rho, rhoFreq) ->
              L.foldl'
                (\(!vec2) ((!phi, !phiFreq), inputVec) ->
                   VS.zipWith
                     (+)
                     vec2
                     (VS.map
                        (* (coefficients R.!
                            (Z :. (0 :: Int) :. theta :. rho :. phi))) .
                      VS.zipWith
                        (*)
                        (getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq 0) $
                      inputVec))
                vec1 .
              L.zip idxTheta $
              dftVecs)
           initVec .
         L.zip [0 ..] $
         rFreqs) $
    idxTheta        
  
  
{-# INLINE convolveG' #-}
convolveG' ::
     Field
  -> DFTPlan
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> VS.Vector (Complex Double)
  -> DFTArray
  -> IO DFTArray
convolveG' !field !plan !coefficients !harmonicsArray !gaussianFilter !arr@(DFTArray rows cols thetaFreqs rFreqs vecs) = do
  dftVecs <-
    L.map (VS.zipWith (*) gaussianFilter) <$>
    dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) vecs
  let !initVec = VS.replicate (VS.length . L.head $ vecs) 0
      idxTheta = L.zip [0 ..] thetaFreqs
  fmap (DFTArray rows cols thetaFreqs rFreqs) .
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) .
    parMap
      rdeepseq
      (\(!theta, !thetaFreq) ->
         L.foldl'
           (\vec1 (rho, rhoFreq) ->
              L.foldl'
                (\(!vec2) ((!phi, !phiFreq), inputVec) ->
                   VS.zipWith
                     (+)
                     vec2
                     (VS.map
                        (* (case field of
                              Source ->
                                coefficients R.!
                                (Z :. (0 :: Int) :. theta :. rho :. phi)
                              Sink ->
                                (coefficients R.!
                                 (Z :. (0 :: Int) :. theta :. rho :. phi)) *
                                (cis $ (-(phiFreq + thetaFreq)) * pi))) .
                      VS.zipWith
                        (*)
                        (getHarmonics harmonicsArray phiFreq rhoFreq thetaFreq 0) $
                      inputVec))
                vec1 .
              L.zip idxTheta $
              dftVecs)
           initVec .
         L.zip [0 ..] $
         rFreqs) $
    idxTheta        

-- {-# INLINE convolve'' #-}
-- convolve'' ::
--      Field
--   -> DFTPlan
--   -> R.Array U DIM4 (Complex Double)
--   -> IA.Array (Int, Int) (VS.Vector (Complex Double))
--   -> DFTArray
--   -> IO DFTArray
-- convolve'' !field !plan !coefficients !harmonicsArray !arr@(DFTArray rows cols phiFreqs rhoFreqs vecs) = do
--   let !initVec = VS.replicate (VS.length . L.head $ vecs) 0
--       idx = L.zip [0 ..] phiFreqs
--   return .
--     (DFTArray rows cols phiFreqs rhoFreqs) .
--     parMap
--       rdeepseq
--       (\(theta, thetaFreq) ->
--          L.foldl'
--            (\vec1 (rho, rhoFreq) ->
--               L.foldl'
--                 (\vec2 (phi, phiFreq) ->
--                    VS.zipWith (+) vec2 .
--                    VS.map
--                      (* (case field of
--                            Source ->
--                              coefficients R.!
--                              (Z :. (0 :: Int) :. theta :. rho :. phi)
--                            Sink ->
--                              (coefficients R.!
--                               (Z :. (0 :: Int) :. theta :. rho :. phi)) *
--                              (cis $ (-(thetaFreq + phiFreq)) * pi))) $
--                    (getHarmonics harmonicsArray phiFreq rhoFreq 0 0))
--                 vec1
--                 idx)
--            initVec .
--          L.zip [0 ..] $
--          rhoFreqs) $
--     idx

-- {-# INLINE dftHarmonicsArrayGPU #-}
-- dftHarmonicsArrayGPU ::
--      DFTPlan
--   -> Int
--   -> Double
--   -> Int
--   -> Double
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> [Double]
--   -> Double
--   -> Double
--   -> IO (Acc (A.Array A.DIM4 (A.Complex Double)))
-- dftHarmonicsArrayGPU !plan !numRows !deltaRow !numCols !deltaCol !phiFreqs !rhoFreqs !thetaFreqs !rFreqs !halfLogPeriod !cutoff = do
--   let !arr =
--         computeHarmonicsArrayGPU
--           numRows
--           deltaRow
--           numCols
--           deltaCol
--           phiFreqs
--           rhoFreqs
--           thetaFreqs
--           rFreqs
--           halfLogPeriod
--           cutoff
--       rangeFunc2 xs ys =
--         (round (L.head xs - L.last ys), round (L.last xs - L.head ys))
--       (!tfLB, !tfUB) = rangeFunc2 phiFreqs thetaFreqs
--       (!rfLB, !rfUB) = rangeFunc2 rhoFreqs rFreqs
--   fmap
--     (A.use .
--      A.fromList
--        (A.Z A.:. (rfUB - rfLB + 1) A.:. (tfUB - tfLB + 1) A.:. numCols A.:.
--         numRows) .
--      VS.toList . VS.concat) .
--     dftExecuteBatchP plan (DFTPlanID DFT1DG [numCols, numRows] [0, 1]) .
--     L.map
--       (VS.convert .
--        toUnboxed .
--        computeUnboxedS . makeFilter2D . fromUnboxed (Z :. numCols :. numRows)) .
--     IA.elems $
    -- arr

-- {-# INLINE convolveGPU #-}
-- convolveGPU ::
--      PTX
--   -> Field
--   -> DFTPlan
--   -> Acc (A.Array A.DIM4 (A.Complex Double))
--   -> Acc (A.Array A.DIM4 (A.Complex Double))
--   -> Acc (A.Array A.DIM4 (A.Complex Double))
--   -> DFTArray
--   -> IO DFTArray
-- convolveGPU !ptx !field !plan !coefficients !coefficientsSink !harmonics !arr@(DFTArray rows cols thetaFreqs rFreqs vecs) = do
--   dftVecs <- dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) vecs
--   let !numRFreq = L.length rFreqs
--       !numThetaFreq = L.length thetaFreqs
--       idx =
--         [ (r, theta)
--         | r <- [0 .. (numRFreq - 1)]
--         , theta <- [0 .. (numThetaFreq - 1)]
--         ]
--   case field of
--     Source ->
--       fmap (DFTArray rows cols thetaFreqs rFreqs) .
--       dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) .
--       L.map
--         (\(r, theta) ->
--            VS.fromList . A.toList $
--            (runNWith ptx (convolveKernel coefficients harmonics))
--              (A.fromList A.Z [theta])
--              (A.fromList A.Z [r])
--              (A.fromList
--                 (A.Z A.:. numRFreq A.:. numThetaFreq A.:. cols A.:. rows) .
--               VS.toList . VS.concat $
--               dftVecs)) $
--       idx
--     -- Sink ->
--     --   fmap (DFTArray rows cols thetaFreqs rFreqs) .
--     --   dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) .
--     --   L.map
--     --     (\(r, theta) ->
--     --        VS.fromList .
--     --        A.toList .
--     --        (runNWith ptx (convolveKernel coefficientsSink harmonics))
--     --          (A.fromList (Z :. (1 :: Int)) theta)
--     --          (A.fromList (Z :. (1 :: Int)) r)
--     --          (A.fromList (Z :. numRFreq :. numThetaFreq :. cols :. rows) .
--     --           VS.toList . VS.concat $
--     --           dftVecs)) $
--     --   idx


{-# INLINE convolveSingle #-}
convolveSingle ::
     Field
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Int
  -> Int
  -> R.Array U DIM2 (Complex Double)
  -> VU.Vector (Complex Double)
convolveSingle !field !coefficients !harmonicsArray !thetaFreqs !rFreqs !x !y !input =
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
      rowCenter = div rows 2
      colCenter = div cols 2
      (Z :. _ :. _ :. numRFreq :. numThetaFreq) = extent coefficients
      harmonics =
        R.traverse2 thetaFreqs rFreqs (\_ _ -> extent coefficients) $ \fThetaFreq fRFreq (Z :. r :. theta :. rho :. phi) ->
          (harmonicsArray IA.!
           ( round (fRFreq (Z :. rho) - fRFreq (Z :. r))
           , round $ (fThetaFreq (Z :. phi) - (fThetaFreq (Z :. theta))))) R.!
          (Z :. (x + colCenter) :. (y + rowCenter))
      product =
        R.traverse2 (R.zipWith (*) coefficients harmonics) input const $ \fCoefficient fInput idx@(Z :. _ :. _ :. rho :. phi) ->
          fCoefficient idx * fInput (Z :. rho :. phi)
      arr =
        case field of
          Source -> product
          Sink ->
            error
              "convolveSingle: You should not use it to compute a sink field."
  in if x == 0 && y == 0
       then VU.replicate (numRFreq * numThetaFreq) 0
       else toUnboxed . sumS . sumS $ arr

-- {-# INLINE convolveSparse #-}
convolveSparse ::
     Field
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> [(Int, Int)]
  -> [R.Array U DIM2 (Complex Double)]
  -> [R.Array U DIM2 (Complex Double)]
convolveSparse !field !coefficients !harmonicsArray !thetaFreqs !rFreqs !cutoff !xs !inputs =
  L.map (fromUnboxed (extent . L.head $ inputs)) .
  parMap
    rdeepseq
    (\(x, y) ->
       L.foldl1' (VU.zipWith (+)) .
       L.zipWith
         (\(x', y') input ->
            let !x'' = x - x'
                !y'' = y - y'
                (!a, !b) =
                  if (fromIntegral $ x'' ^ 2 + y'' ^ 2) <= cutoff ^ 2
                    then (x'', y'')
                    else (0, 0)
            in convolveSingle
                 field
                 coefficients
                 harmonicsArray
                 thetaFreqs
                 rFreqs
                 a
                 b
                 input)
         xs $
       inputs) $
  xs


{-# INLINE convolveSingle' #-}
convolveSingle' ::
     Field
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> Double
  -> R.Array U DIM1 (Complex Double)
  -> VU.Vector (Complex Double)
convolveSingle' !field !coefficients !harmonicsArray !phiFreqs !rhoFreqs !x !y !input =
  let (Z :. cols :. rows) = extent (harmonicsArray IA.! (0, 0))
      rowCenter = div rows 2
      colCenter = div cols 2
      (Z :. _ :. _ :. numRhoFreq :. numPhiFreq) = extent coefficients
      -- harmonics =
      --   R.traverse2 phiFreqs rhoFreqs (\_ _ -> extent coefficients) $ \fPhiFreq fRhoFreq (Z :. _ :. theta :. rho :. phi) ->
      --     (harmonicsArray IA.!
      --      ( round (fRhoFreq (Z :. rho))
      --      , round (fPhiFreq (Z :. phi) - (fPhiFreq (Z :. theta))))) R.!
      --     (Z :. (x + colCenter) :. (y + rowCenter))
      harmonics =
        R.traverse2 phiFreqs rhoFreqs (\_ _ -> extent coefficients) $ \fPhiFreq fRhoFreq (Z :. _ :. theta :. rho :. phi) ->
          let !tf = fPhiFreq (Z :. phi) - fPhiFreq (Z :. theta)
              !rf = fRhoFreq (Z :. rho)
          in (x :+ y) ** (tf :+ 0) *
             ((x ^ 2 + y ^ 2) :+ 0) ** (((-tf - 1) :+ rf) / 2)
      product =
        R.traverse2 (R.zipWith (*) coefficients harmonics) input const $ \fCoefficient fInput idx@(Z :. _ :. _ :. _ :. phi) ->
          fCoefficient idx * fInput (Z :. phi)
      arr =
        case field of
          Source -> product
          Sink ->
            R.traverse2 product phiFreqs const $ \fProd fPhiFreq idx@(Z :. _ :. theta :. _ :. phi) ->
              fProd idx *
              (cis $ (-(fPhiFreq (Z :. phi) - fPhiFreq (Z :. theta))) * pi)
  in if x == 0 && y == 0
       then VU.replicate (numRhoFreq * numPhiFreq) 0
       else toUnboxed . sumS . sumS $ arr


convolveSparse' ::
     Field
  -> R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (R.Array U DIM2 (Complex Double))
  -> R.Array U DIM1 Double
  -> R.Array U DIM1 Double
  -> Double
  -> [(Double, Double)]
  -> [R.Array U DIM1 (Complex Double)]
  -> [R.Array U DIM1 (Complex Double)]
convolveSparse' !field !coefficients !harmonicsArray !phiFreqs !rhoFreqs !cutoff !xs !inputs =
  L.map (fromUnboxed (extent . L.head $ inputs)) .
  parMap
    rdeepseq
    (\(x, y) ->
       L.foldl1' (VU.zipWith (+)) .
       L.zipWith
         (\(x', y') input ->
            let !x'' = x - x'
                !y'' = y - y'
                (!a, !b) =
                  if (x'' ^ 2 + y'' ^ 2) <= cutoff ^ 2
                    then (x'', y'')
                    else (0, 0)
            in convolveSingle'
                 field
                 coefficients
                 harmonicsArray
                 phiFreqs
                 rhoFreqs
                 a
                 b
                 input)
         xs $
       inputs) $
  xs

{-# INLINE gaussianFilter2D #-}
gaussianFilter2D :: DFTPlan -> Int -> Double -> IO (VS.Vector (Complex Double))
gaussianFilter2D plan numPoint stdG =
  dftExecute plan (DFTPlanID DFT1DG [numPoint, numPoint] [0, 1]) .
  VU.convert .
  toUnboxed . computeS . makeFilter2D . fromFunction (Z :. numPoint :. numPoint) $ \(Z :. i :. j) ->
    let x = fromIntegral $ i - div numPoint 2
        y = fromIntegral $ j - div numPoint 2
    in (exp (-(x ^ 2 + y ^ 2) / (2 * stdG ^ 2))) :+ 0


{-# INCLUDE convolveFull' #-}
convolveFull' ::
     R.Array U DIM4 (Complex Double)
  -> IA.Array (Int, Int) (VS.Vector (Complex Double))
  -> DFTArray
  -> DFTArray
convolveFull' !coefficients !harmonicsArray !arr@(DFTArray r2Freq _ thetaFreqs rFreqs vecs) =
  let r2Freqs = [-r2Freq .. r2Freq]
      idxTheta = L.zip [0 ..] thetaFreqs
      !initVec = VS.replicate (VS.length . L.head $ vecs) 0
  in DFTArray r2Freq r2Freq thetaFreqs rFreqs .
     parMap
       rdeepseq
       (\(!theta, !thetaFreq) ->
          L.foldl'
            (\vec1 (rho, rhoFreq) ->
               L.foldl'
                 (\(!vec2) ((!phi, !phiFreq), inputVec) ->
                    VS.zipWith
                      (+)
                      vec2
                      (VS.map
                         (* (coefficients R.!
                             (Z :. (0 :: Int) :. theta :. rho :. phi))) .
                       VS.zipWith
                         (*)
                         (getHarmonics
                            harmonicsArray
                            phiFreq
                            rhoFreq
                            thetaFreq
                            0) $
                       inputVec))
                 vec1 .
               L.zip idxTheta $
               vecs)
            initVec .
          L.zip [0 ..] $
          rFreqs) $
     idxTheta
