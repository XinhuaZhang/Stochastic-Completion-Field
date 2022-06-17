module Image.Edge where

import           Control.Monad       as M
import           Data.Char
import           Data.List           as L
import           Data.Text.Lazy      as TL
import           Data.Text.Lazy.IO   as TL
import           Data.Text.Lazy.Read as TL
import           System.IO
import           Types

{-# INLINE parseInt #-}
parseInt :: Text -> (Int,Text)
parseInt txt =
  case (decimal . TL.dropWhile (not . isDigit) $ txt) of
    Left msg        -> error "parseInt: empty text."
    Right (x, txt1) -> (x, txt1)

{-# INLINE parseDouble #-}
parseDouble :: Text -> (Double,Text)
parseDouble txt =
  case (double . TL.dropWhile (not . isDigit) $ txt) of
    Left msg        -> error "parseDouble: empty text."
    Right (x, txt1) -> (x, txt1)

{-# INLINE parseEdge #-}
parseEdge :: Text -> R2S1RPPoint
parseEdge txt =
  let (x, txt1) = parseInt txt
      (y, txt2) = parseInt txt1
      (theta, _) = parseDouble txt2
   in R2S1RPPoint (x, y, theta, 1)

parseEdgeFile :: FilePath -> IO [R2S1RPPoint]
parseEdgeFile filePath =
  withFile filePath ReadMode $ \h -> do
    numStr <- TL.hGetLine h
    let (num, _) = parseInt numStr
    M.replicateM num $ do
      str <- TL.hGetLine h
      return . parseEdge $ str
      
  
{-# INLINE parseEdge' #-}
parseEdge' :: Text -> (Double, Double)
parseEdge' txt =
  let (x, txt1) = parseDouble txt
      (y, txt2) = parseDouble txt1
      -- (theta, _) = parseDouble txt2
   in (x, y)

parseEdgeFile' :: FilePath -> IO [(Double,Double)]
parseEdgeFile' filePath =
  withFile filePath ReadMode $ \h -> do
    numStr <- TL.hGetLine h
    let (num, _) = parseInt numStr
    M.replicateM num $ do
      str <- TL.hGetLine h
      return . parseEdge' $ str

{-# INLINE parseEdge'' #-}
parseEdge'' :: Text -> Double
parseEdge'' txt =
  let (x, txt1) = parseInt txt
   in fromIntegral x

parseEdgeFile'' :: FilePath -> IO [Double]
parseEdgeFile'' filePath =
  withFile filePath ReadMode $ \h -> do
    numStr <- TL.hGetLine h
    let (num, _) = parseInt numStr
    M.replicateM num $ do
      str <- TL.hGetLine h
      return . parseEdge'' $ str
