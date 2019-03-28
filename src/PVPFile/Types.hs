module PVPFile.Types where

import           Control.DeepSeq
import           Data.Array.Repa     as R
import           Data.Vector.Unboxed as VU

data PVPHeader = PVPHeader
  { headerSize   :: Int
  , numParams    :: Int
  , fileType     :: Int
  , nx           :: Int
  , ny           :: Int
  , nf           :: Int
  , numRecords   :: Int
  , recordSize   :: Int
  , dataSize     :: Int
  , dataType     :: Int
  , nxProcs      :: Int
  , nyProcs      :: Int
  , nxGlobal     :: Int
  , nyGlobal     :: Int
  , kx           :: Int
  , ky           :: Int
  , nb           :: Int
  , nBands       :: Int
  , time         :: Double
  , weightHeader :: PVPWeightHeader
  } deriving (Show)

data PVPWeightHeader = PVPWeightHeader
  { nxp        :: Int
  , nyp        :: Int
  , nfp        :: Int
  , wMin       :: Double
  , wMax       :: Double
  , numPatches :: Int
  } deriving (Show)

data PVPFileType
  = PVP_FILE
  | PVP_ACT_FILE
  | PVP_WGT_FILE
  | PVP_NONSPIKING_ACT_FILE
  | PVP_KERNEL_FILE
  | PVP_ACT_SPARSEVALUES_FILE
  deriving (Show, Eq)

data PVPDataType
  = PV_BYTE
  | PV_INT
  | PV_FLOAT
  | PV_SPARSEVALUES
  deriving (Show)

data PVPDimension = PVPDimension
  { outputNx :: !Int
  , outputNy :: !Int
  , outputNf :: !Int
  } deriving (Show)

instance NFData PVPDimension where
  rnf (PVPDimension a b c) = a `seq` b `seq` c `seq` ()

data PVPOutputData
  = PVP_OUTPUT_ACT !PVPDimension
                   ![Int]
  | PVP_OUTPUT_NONSPIKING_ACT !PVPDimension
                              !(VU.Vector Double)
  | PVP_OUTPUT_ACT_SPARSEVALUES !PVPDimension
                                !(VU.Vector (Int, Double))
  | PVP_OUTPUT_KERNEL !(Array U DIM4 Double)

instance Show PVPOutputData where
  show (PVP_OUTPUT_ACT _ _) = "PVP_OUTPUT_ACT"
  show (PVP_OUTPUT_NONSPIKING_ACT _ _) = "PVP_OUTPUT_NONSPIKING_ACT"
  show (PVP_OUTPUT_ACT_SPARSEVALUES _ _) = "PVP_OUTPUT_ACT_SPARSEVALUES"
  show (PVP_OUTPUT_KERNEL _) = "PVP_OUTPUT_KERNEL"


instance NFData PVPOutputData where
  rnf x =
    case x of
      PVP_OUTPUT_ACT d xs              -> d `seq` xs `seq` ()
      PVP_OUTPUT_NONSPIKING_ACT d xs   -> d `seq` xs `seq` ()
      PVP_OUTPUT_ACT_SPARSEVALUES d xs -> d `seq` xs `seq` ()
      PVP_OUTPUT_KERNEL arr            -> deepSeqArray arr ()

data PVPFrameData
  = FRAME_ACT !Int
  | FRAME_NONSPIKING_ACT !Double
  | FRAME_ACT_SPARSEVALUES !(Int, Double)
  | FRAME_KERNEL ![Double]
  deriving (Show)
