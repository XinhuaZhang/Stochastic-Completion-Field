module Utils.Time where

import           Data.Time
import           Text.Printf

{-# INLINE printCurrentTime #-}
printCurrentTime :: String -> IO ()
printCurrentTime s = do
  time <- getZonedTime
  printf
    "%s %s\n"
    (take 8 . show . localTimeOfDay . zonedTimeToLocalTime $ time)
    s
