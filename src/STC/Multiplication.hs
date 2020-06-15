{-# LANGUAGE Strict #-}
module STC.Multiplication where

import           Control.Concurrent.Async
import           Control.Monad            as M
import           Data.List                as L
import           Data.Vector.Storable     as VS
import           DFT.Plan
import           STC.DFTArray
import           Utils.List
import           Utils.Parallel

{-# INLINE multiplyPinwheelBasis #-}
multiplyPinwheelBasis :: DFTPlan -> DFTArray -> DFTArray -> IO DFTArray
multiplyPinwheelBasis plan (DFTArray rows cols thetaFreqs rFreqs dftVecs1) (DFTArray _ _ _ _ vecs2) = do
  dftVecs2 <- dftExecuteBatchP plan (DFTPlanID DFT1DG [cols, rows] [0, 1]) vecs2
  let dftVecs3 = L.zipWith (VS.zipWith (*)) dftVecs1 dftVecs2
  DFTArray rows cols thetaFreqs rFreqs <$>
    dftExecuteBatchP plan (DFTPlanID IDFT1DG [cols, rows] [0, 1]) dftVecs3

{-# INLINE multiplyPinwheelBasisBatch #-}
multiplyPinwheelBasisBatch ::
     DFTPlan -> Int -> DFTArray -> DFTArray -> IO DFTArray
multiplyPinwheelBasisBatch plan numBatch (DFTArray rows cols thetaFreqs rFreqs dftVecs1) (DFTArray _ _ _ _ vecs2) = do
  let forwardPlanID = DFTPlanID DFT1DG [cols, rows] [0, 1]
      backwardPlanID = DFTPlanID IDFT1DG [cols, rows] [0, 1]
      forwardPlan = getDFTPlan plan forwardPlanID
      backwardPlan = getDFTPlan plan backwardPlanID
  fmap (DFTArray rows cols thetaFreqs rFreqs . L.concat) .
    mapConcurrently
      (M.mapM
         (\(dftVec1, vec2) -> do
            dftVec2 <- dftExecuteWithPlan forwardPlanID forwardPlan vec2
            dftExecuteWithPlan backwardPlanID backwardPlan $
              VS.zipWith (*) dftVec1 dftVec2)) .
    divideListN numBatch $
    L.zip dftVecs1 vecs2
