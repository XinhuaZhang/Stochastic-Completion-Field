module STC.Point where

data Point =
  Point {-# UNPACK #-}!Int    -- x
        {-# UNPACK #-}!Int    -- y
        {-# UNPACK #-}!Double --theta
        {-# UNPACK #-}!Double --scale
  deriving (Show, Read)
