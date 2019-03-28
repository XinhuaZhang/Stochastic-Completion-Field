module PVPFile.IO
  ( module PVPFile.Types
  , PVPHeader(..)
  , PVPWeightHeader(..)
  , PVPOutputData(..)
  , PVPDimension(..)
  , readPVPHeader
  , pvpFileSource
  , pvpOutputData2Array
  , hPutPVPHeader
  , writePVPFileConduit
  ) where

import           Control.Monad                 as M
import           Control.Monad.IO.Class
import           Control.Monad.Trans.Resource
import           Data.Array.Repa               as R
import           Data.Binary.Get
import           Data.Binary.Put
import qualified Data.ByteString               as BS
import qualified Data.ByteString.Lazy          as BL
import qualified Data.ByteString.Lazy.Internal as BL
import           Data.Conduit                  as C
import           Data.Conduit.Binary           as CB
import           Data.List                     as L
import           Data.Vector.Unboxed           as VU
import           GHC.Float
import           PVPFile.Types
import           Prelude                       as P
import           System.IO

getPVPFileType :: PVPHeader -> PVPFileType
getPVPFileType header =
  case fileType header of
    1 -> PVP_FILE
    2 -> PVP_ACT_FILE
    3 -> PVP_WGT_FILE
    4 -> PVP_NONSPIKING_ACT_FILE
    5 -> PVP_KERNEL_FILE
    6 -> PVP_ACT_SPARSEVALUES_FILE
    _ -> error "Wrong PVPFileType."

getPVPDataType :: PVPHeader -> PVPDataType
getPVPDataType header =
  case dataType header of
    1 -> PV_BYTE
    2 -> PV_INT
    3 -> PV_FLOAT
    4 -> PV_SPARSEVALUES
    _ -> error "Wrong PVPDataType."

-- Header functions
getHeaderParam :: Get PVPHeader
getHeaderParam = do
  headerSize' <- getWord32le
  numParams' <- getWord32le
  fileType' <- getWord32le
  nx' <- getWord32le
  ny' <- getWord32le
  nf' <- getWord32le
  numRecords' <- getWord32le
  recordSize' <- getWord32le
  dataSize' <- getWord32le
  dataType' <- getWord32le
  nxProcs' <- getWord32le
  nyProcs' <- getWord32le
  nxGlobal' <- getWord32le
  nyGlobal' <- getWord32le
  kx' <- getWord32le
  ky' <- getWord32le
  nb' <- getWord32le
  nBands' <- getWord32le
  time' <- getDoublele
  wHeader' <-
    if fileType' == 3 || fileType' == 5
      then getWeightHeaderParams
      else return $ PVPWeightHeader 0 0 0 0 0 0
  return $
    PVPHeader
      (fromIntegral headerSize')
      (fromIntegral numParams')
      (fromIntegral fileType')
      (fromIntegral nx')
      (fromIntegral ny')
      (fromIntegral nf')
      (fromIntegral numRecords')
      (fromIntegral recordSize')
      (fromIntegral dataSize')
      (fromIntegral dataType')
      (fromIntegral nxProcs')
      (fromIntegral nyProcs')
      (fromIntegral nxGlobal')
      (fromIntegral nyGlobal')
      (fromIntegral kx')
      (fromIntegral ky')
      (fromIntegral nb')
      (fromIntegral nBands')
      time'
      wHeader'

getWeightHeaderParams :: Get PVPWeightHeader
getWeightHeaderParams = do
  nxp' <- getInt32le
  nyp' <- getInt32le
  nfp' <- getInt32le
  wMin' <- getFloatle
  wMax' <- getFloatle
  numPatches' <- getInt32le
  return $!
    PVPWeightHeader
      (fromIntegral nxp')
      (fromIntegral nyp')
      (fromIntegral nfp')
      (float2Double wMin')
      (float2Double wMax')
      (fromIntegral numPatches')

getPVPHeader :: Handle -> IO PVPHeader
getPVPHeader h = do
  bs <- BL.hGet h 4
  let headerSize' = fromIntegral $ runGet getInt32le bs
  bs' <- BL.hGet h (headerSize' - 4)
  return $ runGet getHeaderParam (BL.append bs bs')

readPVPHeader :: FilePath -> IO PVPHeader
readPVPHeader filePath = do
  h <- openBinaryFile filePath ReadMode
  header <- getPVPHeader h
  hClose h
  return header

getPVPDimension :: PVPHeader -> PVPDimension
getPVPDimension header = PVPDimension (nx header) (ny header) (nf header)

-- Frame functions
getByteStringData :: Handle -> Int -> PVPDataType -> IO BL.ByteString
getByteStringData handle n dataType' =
  case dataType' of
    PV_BYTE         -> BL.hGet handle n
    PV_INT          -> BL.hGet handle (4 * n)
    PV_FLOAT        -> BL.hGet handle (4 * n)
    PV_SPARSEVALUES -> BL.hGet handle (4 * n)

getPVPFrameData :: PVPHeader -> Get PVPFrameData
getPVPFrameData header =
  case getPVPFileType header of
    PVP_FILE -> undefined
    PVP_ACT_FILE -> do
      ind <- getWord32le
      return $ FRAME_ACT (fromIntegral ind)
    PVP_WGT_FILE -> undefined
    PVP_NONSPIKING_ACT_FILE -> do
      val <- getFloatle
      return $ FRAME_NONSPIKING_ACT (float2Double val)
    PVP_KERNEL_FILE -> do
      _ <- getInt16le
      _ <- getInt16le
      _ <- getInt32le
      xs <- M.replicateM (nxp' * nyp' * nfp') getFloatle
      return $ FRAME_KERNEL . P.map float2Double $ xs
    PVP_ACT_SPARSEVALUES_FILE -> do
      ind <- getWord32le
      v <- getFloatle
      return $ FRAME_ACT_SPARSEVALUES (fromIntegral ind, float2Double v)
  where
    (PVPWeightHeader nxp' nyp' nfp' _ _ _numPatches') = weightHeader header

incrementalGetPVPFrameData :: PVPHeader -> BL.ByteString -> [PVPFrameData]
incrementalGetPVPFrameData header = go decoder
  where
    decoder = runGetIncremental (getPVPFrameData header)
    go :: Decoder PVPFrameData -> BL.ByteString -> [PVPFrameData]
    go (Done leftover' _consumed x) input
      | BS.null leftover' && BL.null input = [x]
      | otherwise = x : go decoder (BL.chunk leftover' input)
    go (Partial k) input = go (k . takeHeadChunk $ input) (dropHeadChunk input)
    go (Fail _leftover _consumed msg) _input = error msg

takeHeadChunk :: BL.ByteString -> Maybe BS.ByteString
takeHeadChunk lbs =
  case lbs of
    (BL.Chunk bs _) -> Just bs
    _               -> Nothing

dropHeadChunk :: BL.ByteString -> BL.ByteString
dropHeadChunk lbs =
  case lbs of
    (BL.Chunk _ lbs') -> lbs'
    _                 -> BL.Empty

getFrame :: PVPHeader -> Handle -> IO PVPOutputData
getFrame header h =
  case getPVPFileType header of
    PVP_FILE -> undefined
    PVP_ACT_FILE -> do
      bs <- BL.hGet h 12
      let (_time, numActive) =
            runGet
              (do time' <- getDoublele
                  num <- getWord32le
                  return (time', fromIntegral num))
              bs
      bs' <- getByteStringData h numActive (getPVPDataType header)
      return .
        PVP_OUTPUT_ACT (getPVPDimension header) . P.map (\(FRAME_ACT x) -> x) $
        incrementalGetPVPFrameData header bs'
    PVP_WGT_FILE -> undefined
    PVP_NONSPIKING_ACT_FILE -> do
      bs <- BL.hGet h 8
      let _time = runGet getDoublele bs
      bs' <- getByteStringData h (recordSize header) (getPVPDataType header)
      return .
        PVP_OUTPUT_NONSPIKING_ACT (getPVPDimension header) .
        VU.fromList . P.map (\(FRAME_NONSPIKING_ACT x) -> x) $
        incrementalGetPVPFrameData header bs'
    PVP_KERNEL_FILE -> do
      _frameHeader <- getPVPHeader h
      let rs =
            (numPatches . weightHeader $ header) *
            ((nxp . weightHeader $ header) * (nyp . weightHeader $ header) *
             (nfp . weightHeader $ header) *
             (dataSize header) +
             8)
      bs' <- BL.hGet h rs -- (recordSize header)
      let xs =
            P.concat . L.transpose . P.map (\(FRAME_KERNEL x) -> x) $
            incrementalGetPVPFrameData header bs'
      return .
        PVP_OUTPUT_KERNEL .
        fromListUnboxed (Z :. nyp' :. nxp' :. nfp' :. numPatches') $
        xs
    PVP_ACT_SPARSEVALUES_FILE -> do
      bs <- BL.hGet h 12
      let (_time, numActive) =
            runGet
              (do time' <- getDoublele
                  num <- getWord32le
                  return (time', fromIntegral num))
              bs
      if numActive == 0
        then do
          putStrLn "All zeros activities."
          return (PVP_OUTPUT_ACT_SPARSEVALUES (getPVPDimension header) VU.empty)
        else do
          bs' <- getByteStringData h (numActive * 2) (getPVPDataType header)
          return .
            PVP_OUTPUT_ACT_SPARSEVALUES (getPVPDimension header) .
            VU.fromList . P.map (\(FRAME_ACT_SPARSEVALUES x) -> x) $
            incrementalGetPVPFrameData header bs'
  where
    (PVPWeightHeader nxp' nyp' nfp' _ _ numPatches') = weightHeader header

pvpFileSource :: FilePath -> ConduitT () PVPOutputData (ResourceT IO) ()
pvpFileSource filePath = do
  h <- liftIO $ openBinaryFile filePath ReadMode
  header <- liftIO $ getPVPHeader h
  let fileType' = getPVPFileType header
  h1 <-
    if fileType' == PVP_WGT_FILE || fileType' == PVP_KERNEL_FILE
      then do
        liftIO $ hClose h
        liftIO $ openBinaryFile filePath ReadMode
      else return h
  loop (nBands header) h1 header
  where
    loop n handle header' =
      if n > 0
        then do
          frame <- liftIO $ getFrame header' handle
          yield frame
          loop (n - 1) handle header'
        else do
          liftIO $ hClose handle
          return ()

-- layout: ny * nx * nf
pvpOutputData2Array :: PVPOutputData -> Array U DIM3 Double
pvpOutputData2Array (PVP_OUTPUT_NONSPIKING_ACT (PVPDimension nx' ny' nf') vec) =
  fromUnboxed (Z :. ny' :. nx' :. nf') vec
pvpOutputData2Array (PVP_OUTPUT_ACT_SPARSEVALUES (PVPDimension nx' ny' nf') vec) =
  fromUnboxed (Z :. ny' :. nx' :. nf') .
  accumulate (+) (VU.replicate (ny' * nx' * nf') 0) $vec
pvpOutputData2Array _ =
  error "pvpOutputData2Array: pvpOutput format is not supported."

-- Write a PVP file
putPVPHeader :: PVPHeader -> Put
putPVPHeader header =
  case fileType header of
    4 -> do
      putWord32le . fromIntegral . headerSize $ header
      putWord32le . fromIntegral . numParams $ header
      putWord32le . fromIntegral . fileType $ header
      putWord32le . fromIntegral . nx $ header
      putWord32le . fromIntegral . ny $ header
      putWord32le . fromIntegral . nf $ header
      putWord32le . fromIntegral . numRecords $ header
      putWord32le . fromIntegral . recordSize $ header
      putWord32le . fromIntegral . dataSize $ header
      putWord32le . fromIntegral . dataType $ header
      putWord32le . fromIntegral . nxProcs $ header
      putWord32le . fromIntegral . nyProcs $ header
      putWord32le . fromIntegral . nxGlobal $ header
      putWord32le . fromIntegral . nyGlobal $ header
      putWord32le . fromIntegral . kx $ header
      putWord32le . fromIntegral . ky $ header
      putWord32le . fromIntegral . nb $ header
      putWord32le . fromIntegral . nBands $ header
      putDoublele . time $ header
    5 -> do
      putWord32le . fromIntegral . headerSize $ header
      putWord32le . fromIntegral . numParams $ header
      putWord32le . fromIntegral . fileType $ header
      putWord32le . fromIntegral . nx $ header
      putWord32le . fromIntegral . ny $ header
      putWord32le . fromIntegral . nf $ header
      putWord32le . fromIntegral . numRecords $ header
      putWord32le . fromIntegral . recordSize $ header
      putWord32le . fromIntegral . dataSize $ header
      putWord32le . fromIntegral . dataType $ header
      putWord32le . fromIntegral . nxProcs $ header
      putWord32le . fromIntegral . nyProcs $ header
      putWord32le . fromIntegral . nxGlobal $ header
      putWord32le . fromIntegral . nyGlobal $ header
      putWord32le . fromIntegral . kx $ header
      putWord32le . fromIntegral . ky $ header
      putWord32le . fromIntegral . nb $ header
      putWord32le . fromIntegral . nBands $ header
      putDoublele . time $ header
      putWeightHeaderParams . weightHeader $ header
    _ ->
      error $
      "putPVPHear: Doesn't support file type " L.++ show (getPVPFileType header)
      
putWeightHeaderParams :: PVPWeightHeader -> Put
putWeightHeaderParams weightHeader' = do
  putInt32le . fromIntegral . nxp $ weightHeader'
  putInt32le . fromIntegral . nyp $ weightHeader'
  putInt32le . fromIntegral . nfp $ weightHeader'
  putFloatle . double2Float . wMin $ weightHeader'
  putFloatle . double2Float . wMax $ weightHeader'
  putInt32le . fromIntegral . numPatches $ weightHeader'


hPutPVPHeader :: Handle -> PVPHeader -> IO ()
hPutPVPHeader h = BL.hPut h . runPut . putPVPHeader

putPVPFrame :: PVPOutputData -> Put
putPVPFrame (PVP_OUTPUT_NONSPIKING_ACT _ vec) = do
  putDoublele 0 -- time
  VU.mapM_ (putFloatle . double2Float) vec
putPVPFrame (PVP_OUTPUT_KERNEL arr) = do
  let (Z :. nx :. ny :. _ :. np) = extent arr
  M.mapM_
    (\i -> do
       putInt16le . fromIntegral $ nx
       putInt16le . fromIntegral $ ny
       putInt32le 0
       VU.mapM_ (putFloatle . double2Float) . toUnboxed . computeS . R.slice arr $
         (Z :. All :. All :. All :. i))
    [0 .. np - 1]

putPVPFrame x =
  error $ "putPVPFrame: PVP data type " L.++ show x L.++ " is not supported."

writePVPFileConduit :: ConduitT PVPOutputData BS.ByteString (ResourceT IO) ()
writePVPFileConduit = awaitForever (CB.sourceLbs . runPut . putPVPFrame)
