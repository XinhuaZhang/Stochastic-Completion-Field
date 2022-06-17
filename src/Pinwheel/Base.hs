{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Strict #-}
module Pinwheel.Base where

import           Data.Complex

-- {-# INLINE fourierMellin #-}
-- fourierMellin :: (Eq e, RealFloat e) => e -> Int -> Int -> (e, e) -> Complex e
-- fourierMellin sigma angularFreq radialFreq (!x, !y) =
--   if x == 0 && y == 0
--     then 0
--     else (x :+ y) ** (fromIntegral (-angularFreq) :+ 0) *
--          ((x ^ 2 + y ^ 2) :+ 0) **
--          (((fromIntegral angularFreq - 2 + sigma) :+ (-fromIntegral radialFreq)) /
--           2)
          
{-# INLINE fourierMellin #-}
fourierMellin :: (Eq e, RealFloat e) => e -> Int -> Int -> (e, e) -> Complex e
fourierMellin sigma angularFreq radialFreq (!x, !y) =
  if x == 0 && y == 0
    then 0
    else let r = sqrt $ x ^ 2 + y ^ 2
             theta = atan2 y x
         in ((r :+ 0) ** ((sigma - 1) :+ fromIntegral (-radialFreq))) *
            (cis $ fromIntegral (-angularFreq) * theta)
            
{-# INLINE fourierMellinInv #-}
fourierMellinInv :: (Eq e, RealFloat e) => e -> Int -> Int -> (e, e) -> Complex e
fourierMellinInv sigma angularFreq radialFreq (!x, !y) =
  if x == 0 && y == 0
    then 0
    else -- (x :+ y) ** (fromIntegral angularFreq :+ 0) *
         -- ((x ^ 2 + y ^ 2) :+ 0) **
         -- (((-(fromIntegral angularFreq) + sigma) :+ fromIntegral radialFreq) /
         --  2)
                  let r = sqrt $ x^2 + y^2
                      theta = atan2 y x
                  in ((r :+ 0) ** (sigma :+ fromIntegral radialFreq)) * cis (fromIntegral angularFreq * theta)


{-# INLINE fourierMellinInvPeriod #-}
fourierMellinInvPeriod ::
     (Eq e, RealFloat e) => e -> e -> Int -> Int -> (e, e) -> Complex e
fourierMellinInvPeriod sigma periodEnv angularFreq radialFreq (!x, !y) =
  if x == 0 && y == 0
    then 0
    else (x :+ y) ** (fromIntegral angularFreq :+ 0) *
         ((x ^ 2 + y ^ 2) :+ 0) **
         (((-(fromIntegral angularFreq) + sigma) :+ 2 * pi * fromIntegral radialFreq / (log periodEnv)) /
          2)

