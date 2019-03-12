module Resize where

import           Data.Array.Repa    as R
import           Image.IO
import           Image.Transform
import           System.Directory
import           System.Environment
import           System.FilePath
import Data.Vector.Unboxed as VU

main = do
  (inputPath:outputPath:sizeStr:_) <- getArgs
  let size = read sizeStr :: Int
      folderPath = "output/test/Resize"
  (ImageRepa _ img) <- readImageRepa inputPath False
  createDirectoryIfMissing True folderPath
  let arr = resize25D (size, size) (0, 255) img
      avg = R.sumAllS arr / (fromIntegral $ size ^ 2)
  plotImageRepa (folderPath </> outputPath) .
    ImageRepa 8 .
    computeS .
    R.map
      (\x ->
         if x > avg
           then 0
           else 255) $
    arr
