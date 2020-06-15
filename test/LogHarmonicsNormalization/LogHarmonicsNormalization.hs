module LogHarmonicsNormalization where

import Data.Complex
import Data.List as L

main = do
  let delta = 0.00001 :: Double
      maxScale = 16 :: Double
      freq = 5
      xs =
        L.map
          (\x -> cis $ (-pi) * freq * (log x) /  (log maxScale))
          [1 / maxScale,1 / maxScale + delta .. maxScale]
      y =
        (exp (log maxScale :+ (-pi) * freq) -
         exp (((-1) * log maxScale) :+ pi * freq)) / (1 :+ (-pi) * freq / (log maxScale))
  print . polar $ ((delta :+ 0) * (L.sum xs))
  print . polar $ y
