{-# LANGUAGE Strict #-}
module PlotGabor where

import           Data.Array.Repa    as R
import           Image.IO
import           System.Directory
import           System.Environment
import           System.FilePath
import           Text.Printf
import           Utils.Distribution

{-# INLINE rotate #-}
rotate :: Double -> Double -> Double -> (Double,Double)
rotate deltaTheta x y =
  let r = sqrt $ x ^ 2 + y ^ 2
      theta = deltaTheta + atan2 y x
   in (r * cos theta, r * sin theta)

gaborFunction2D ::
     Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Double
  -> Int
  -> IO (R.Array U DIM2 Double)
gaborFunction2D freq stdX stdY delta orientation scale centerX centerY numPoints = do
  let centerR2 = div numPoints 2
      arr =
        R.fromFunction (Z :. numPoints :. numPoints) $ \(Z :. i :. j) ->
          let x = fromIntegral (i - centerR2) * delta / scale
              y = fromIntegral (j - centerR2) * delta / scale
              (x1, y1) =
                rotate orientation (x - centerX / scale) (y - centerY / scale)
              grating = sin (freq * y1)
              envelope = gaussian1D x1 stdX * gaussian1D y1 stdY
           in grating * envelope
  computeUnboxedP arr


main = do
  args@(freqStr:stdXStr:stdYStr:deltaStr:orientationStr:scaleStr:centerXStr:centerYStr:numPointsStr:_) <-
    getArgs
  let freq = read freqStr :: Double
      stdX = read stdXStr :: Double
      stdY = read stdYStr :: Double
      delta = read deltaStr :: Double
      orientation = read orientationStr :: Double
      scale = read scaleStr :: Double
      centerX = read centerXStr :: Double
      centerY = read centerYStr :: Double
      numPoints = read numPointsStr :: Int
      folderPath = "output/test/PlotGabor"
  createDirectoryIfMissing True folderPath
  gaborArr <-
    gaborFunction2D
      freq
      stdX
      stdY
      delta
      (orientation / 180 * pi)
      scale
      centerX
      centerY
      numPoints
  plotImageRepa
    (folderPath </>
     printf
       "Gabor_%d_%d_%d_%.2f.png"
       (round centerX :: Int)
       (round centerY :: Int)
       (round orientation :: Int)
       scale) .
    ImageRepa 8 . computeS . extend (Z :. (1 :: Int) :. All :. All) $
    gaborArr
