module STC.Point where

data Point =
  Point {-# UNPACK #-}!Double    -- x
        {-# UNPACK #-}!Double    -- y
        {-# UNPACK #-}!Double --theta
        {-# UNPACK #-}!Double --scale
  deriving (Show, Read)
  
