module Resize where

import           Image.IO
import           Image.Transform
import           System.Directory
import           System.Environment
import System.FilePath

main = do
  (inputPath:outputPath:sizeStr:_) <- getArgs
  let size = read sizeStr :: Int
      folderPath = "output/test/Resize"
  (ImageRepa _ img) <- readImageRepa inputPath False
  createDirectoryIfMissing True folderPath
  plotImageRepa (folderPath </> outputPath) .
    ImageRepa 8 . resize25D (size, size) (0, 255) $
    img
