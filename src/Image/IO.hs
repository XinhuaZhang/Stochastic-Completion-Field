module Image.IO
  ( ImageRepa(..)
  , readImagePathList
  , imagePathSource
  , readImageRepa
  , readImageConduit
  , plotImageRepa
  , complexImageToColorImage
  , plotImageRepaComplex
  , plotImageRepaComplexGray
  ) where

import           Codec.Picture                
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Trans.Resource
import           Data.Array.Repa              as R
import           Data.Complex
import           Data.Conduit                 as C
import           Data.Conduit.List            as CL
import           Data.List                    as L
import           Data.Vector.Unboxed          as VU
import           Data.Word
import           GHC.Float

-- x = cols = width
-- y = rows = height
data ImageRepa a = ImageRepa
  { imageDepth   :: !Int
  , imageContent :: !(R.Array U DIM3 a) -- nf x cols x rows
  }

readImagePathList :: FilePath -> IO [String]
readImagePathList = fmap lines . readFile

imagePathSource :: FilePath -> ConduitT () FilePath (ResourceT IO) ()
imagePathSource filePath = do
  pathList <- liftIO $ readImagePathList filePath
  sourceList pathList

readImageRepa :: FilePath -> Bool -> IO (ImageRepa Double)
readImageRepa filePath isColor = do
  buffer <- liftIO $ readImage filePath
  case buffer of
    Left msg -> error msg
    Right dImg ->
      return $
      if isColor
        then case dImg of
               ImageY8 img ->
                 ImageRepa 8 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      fromIntegral $ pixelAt img x y :: Double)
               ImageY16 img ->
                 ImageRepa 16 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      fromIntegral $ pixelAt img x y :: Double)
               ImageYF img ->
                 ImageRepa 64 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      float2Double $ pixelAt img x y :: Double)
               ImageRGB8 img ->
                 ImageRepa 8 . computeS $
                 fromFunction
                   (Z :. (3 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. k :. x :. y) ->
                      let (PixelRGB8 r g b) = pixelAt img x y
                       in case k of
                            0 -> fromIntegral r
                            1 -> fromIntegral g
                            2 -> fromIntegral b
                            _ -> error "readImageConduit: dimension error.")
               ImageRGB16 img ->
                 ImageRepa 16 . computeS $
                 fromFunction
                   (Z :. (3 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. k :. x :. y) ->
                      let (PixelRGB16 r g b) = pixelAt img x y
                       in case k of
                            0 -> fromIntegral r
                            1 -> fromIntegral g
                            2 -> fromIntegral b
                            _ -> error "readImageConduit: dimension error.")
               ImageRGBF img ->
                 ImageRepa 64 . computeS $
                 fromFunction
                   (Z :. (3 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. k :. x :. y) ->
                      let (PixelRGBF r g b) = pixelAt img x y
                       in case k of
                            0 -> float2Double r
                            1 -> float2Double g
                            2 -> float2Double b
                            _ -> error "readImageConduit: dimension error.")
               img ->
                 let rgbImg = convertRGB8 img
                  in ImageRepa 8 . computeS $
                     fromFunction
                       (Z :. (3 :: Int) :. imageWidth rgbImg :.
                        imageHeight rgbImg)
                       (\(Z :. k :. x :. y) ->
                          let (PixelRGB8 r g b) = pixelAt rgbImg x y
                           in case k of
                                0 -> fromIntegral r
                                1 -> fromIntegral g
                                2 -> fromIntegral b
                                _ -> error "readImageConduit: dimension error.")
        else case dImg of
               ImageY8 img ->
                 ImageRepa 8 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      fromIntegral $ pixelAt img x y :: Double)
               ImageY16 img ->
                 ImageRepa 16 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      fromIntegral $ pixelAt img x y :: Double)
               ImageYF img ->
                 ImageRepa 64 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      float2Double $ pixelAt img x y :: Double)
               ImageRGB8 img ->
                 ImageRepa 8 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      let (PixelRGB8 r g b) = pixelAt img x y
                       in rgb2Gray
                            (fromIntegral (maxBound :: Word8))
                            (fromIntegral r)
                            (fromIntegral g)
                            (fromIntegral b))
               ImageRGB16 img ->
                 ImageRepa 16 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      let (PixelRGB16 r g b) = pixelAt img x y
                       in rgb2Gray
                            (fromIntegral (maxBound :: Word16))
                            (fromIntegral r)
                            (fromIntegral g)
                            (fromIntegral b))
               ImageRGBF img ->
                 ImageRepa 64 . computeS $
                 fromFunction
                   (Z :. (1 :: Int) :. imageWidth img :. imageHeight img)
                   (\(Z :. _ :. x :. y) ->
                      let (PixelRGBF r g b) = pixelAt img x y
                       in rgb2Gray
                            1
                            (float2Double r)
                            (float2Double g)
                            (float2Double b))
               img ->
                 let rgbImg = convertRGB8 img
                  in ImageRepa 8 . computeS $
                     fromFunction
                       (Z :. (1 :: Int) :. imageWidth rgbImg :.
                        imageHeight rgbImg)
                       (\(Z :. _ :. x :. y) ->
                          let (PixelRGB8 r g b) = pixelAt rgbImg x y
                           in rgb2Gray
                                (fromIntegral (maxBound :: Word8))
                                (fromIntegral r)
                                (fromIntegral g)
                                (fromIntegral b))

readImageConduit ::
     Bool -> ConduitT FilePath (ImageRepa Double) (ResourceT IO) ()
readImageConduit isColor =
  awaitForever
    (\filePath -> do
       img <- liftIO $ readImageRepa filePath isColor
       yield img)

{-# INLINE rgb2Gray #-}
rgb2Gray :: Double -> Double -> Double -> Double -> Double
rgb2Gray bound r g b
  | yLinear <= 0.0031308 = 12.92 * yLinear * bound
  | otherwise = (1.055 * (yLinear ** (1 / 2.4)) - 0.055) * bound
  where
    yLinear =
      0.2126 * func bound r + 0.7152 * func bound g + 0.0722 * func bound b

{-# INLINE func #-}
func :: Double -> Double -> Double
func bound x
  | y < 0.04045 = y / 12.92
  | otherwise = ((y + 0.055) / 1.055) ** 2.4
  where
    y = x / bound

{-# INLINE normalizeImageRepa #-}
normalizeImageRepa :: (RealFrac e, Eq e, Unbox e, Ord e) => ImageRepa e -> ImageRepa e
normalizeImageRepa i@(ImageRepa depth img)
  | maxV == minV =
    ImageRepa depth . computeS . R.map (\x -> (fromIntegral $ 2 ^ depth - 1)) $
    img
  | otherwise =
    ImageRepa depth .
    computeS .
    R.map (\x -> (x - minV) / (maxV - minV) * (fromIntegral $ 2 ^ depth - 1)) $
    img
  where
    maxV = VU.maximum . toUnboxed $ img
    minV = VU.minimum . toUnboxed $ img

{-# INLINE plotImageRepa #-}
plotImageRepa :: (Unbox e, Ord e, RealFrac e) => FilePath -> ImageRepa e -> IO ()
plotImageRepa filePath img@(ImageRepa depth x) = do
  let Z :. nfp' :. nxp' :. nyp' = extent x
      normalizedImg = imageContent . normalizeImageRepa $ img
      w =
        case nfp' of
          1 ->
            -- ImageY8 $
            -- generateImage
            --   (\x y ->
            --      let v =
            --            fromIntegral . round $
            --            normalizedImg R.! (Z :. 0 :. x :. y)
            --       in v)
            --   nxp'
            --   nyp'
            ImageRGB8 $
            generateImage
              (\x y ->
                 let r =
                       fromIntegral . round $
                       normalizedImg R.! (Z :. 0 :. x :. y)
                     g =
                       fromIntegral . round $
                       normalizedImg R.! (Z :. 0 :. x :. y)
                     b =
                       fromIntegral . round $
                       normalizedImg R.! (Z :. 0 :. x :. y)
                  in PixelRGB8 r g b)
              nxp'
              nyp'
          3 ->
            ImageRGB8 $
            generateImage
              (\x y ->
                 let r =
                       fromIntegral . round $
                       normalizedImg R.! (Z :. 0 :. x :. y)
                     g =
                       fromIntegral . round $
                       normalizedImg R.! (Z :. 1 :. x :. y)
                     b =
                       fromIntegral . round $
                       normalizedImg R.! (Z :. 2 :. x :. y)
                  in PixelRGB8 r g b)
              nxp'
              nyp'
          _ ->
            error $
            "ImageRepa is neither a gray image nor a color image. There are " L.++
            show nfp' L.++
            " channels."
  savePngImage filePath w

-- Plot Complex Image
{-# INLINE complexScale #-}
complexScale :: (Ord e, Floating e, Unbox e) => ImageRepa (Complex e) -> e
complexScale (ImageRepa _ arr) =
  let (rVec, iVec) = VU.unzip . VU.map (\(a :+ b) -> (a, b)) . R.toUnboxed $ arr
      rMax = VU.maximum rVec
      rMin = VU.minimum rVec
      iMax = VU.maximum iVec
      iMin = VU.minimum iVec
   in 2 / ((max rMax iMax) - (min rMin iMin))

{-# INLINE complexImageToColorImage #-}
complexImageToColorImage ::
     (Floating e, Ord e, Unbox e) => ImageRepa (Complex e) -> ImageRepa e
complexImageToColorImage img@(ImageRepa depth arr) =
  let scale = complexScale img
      (Z :. _ :. cols :. rows) = extent arr
      (r, g, b) =
        VU.unzip3 .
        VU.map
          (\(r :+ i) ->
             let [x, y] = L.map (* scale) [r, i]
                 radius = sqrt (x * x + y * y)
                 a = 0.40824829046386301636 * x
                 b = 0.70710678118654752440 * y
                 d = 1 / (1 + radius * radius)
                 d' = 0.5 - radius * d
                 r' = 0.5 + (0.81649658092772603273 * x * d)
                 b' = 0.5 - (d * (a - b))
                 g' = 0.5 - (d * (a + b))
                 [red, grn, blu] = L.map (+ d') [r', g', b']
                 [red', grn', blu'] = L.map (flip (-) d') [r', g', b']
              in if radius < 1
                   then (red', grn', blu')
                   else (red, grn, blu)) .
        toUnboxed $
        arr
   in ImageRepa 8 . fromUnboxed (Z :. (3 :: Int) :. cols :. rows) . VU.concat $
      [r, g, b]
      
{-# INLINE normalizeComplexImage #-}
normalizeComplexImage ::
     ImageRepa (Complex Double) -> ImageRepa (Complex Double)
normalizeComplexImage (ImageRepa depth img) =
  let maxV = VU.maximum magVec
      minV = VU.minimum magVec
      maxVal = 2 ^ depth
      minVal = 0
      magVec = VU.map magnitude . toUnboxed $ img
   in if maxV == minV
        then (ImageRepa depth img)
        else ImageRepa depth .
             fromUnboxed (extent img) .
             VU.zipWith (\a b -> mkPolar b (phase a)) (toUnboxed img) .
             VU.map
               (\x -> (x - minV) / (maxV - minV) * (maxVal - minVal) + minVal) $
             magVec
             

{-# INLINE plotImageRepaComplex #-}
plotImageRepaComplex ::
     (Unbox e, Floating e, RealFrac e)
  => FilePath
  -> ImageRepa (Complex e)
  -> IO ()
plotImageRepaComplex filePath img =
  plotImageRepa filePath . complexImageToColorImage -- . normalizeComplexImage
  $
  img
  
{-# INLINE plotImageRepaComplexGray #-}
plotImageRepaComplexGray :: FilePath -> ImageRepa (Complex Double) -> IO ()
plotImageRepaComplexGray filePath img =
  let (ImageRepa depth arr) = complexImageToColorImage img
      (Z :. chs :. cols :. rows) = extent arr
  in plotImageRepa filePath .
     ImageRepa depth .
     computeS .
     extend (Z :. (1 :: Int) :. All :. All) .
     R.sumS .
     R.backpermute
       (Z :. cols :. rows :. chs)
       (\(Z :. c :. r :. k) -> (Z :. k :. c :. r)) $
     arr
