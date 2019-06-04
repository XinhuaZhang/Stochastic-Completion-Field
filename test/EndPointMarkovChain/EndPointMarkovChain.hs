module EndPointMarkovChain where

import           Control.Monad           as M
import           Data.Array.Repa         as R
import           Data.List               as L
import           Graphics.Gnuplot.Simple
import           Numeric.LinearAlgebra   as NL
import           System.Directory
import           System.Environment
import           System.FilePath         ((</>))
import           Text.Printf

{-# INLINE makeMatrix #-}
makeMatrix :: Int -> Double -> R.Array D DIM2 Double
makeMatrix n p =
  fromFunction (Z :. (2 * n) :. (2 * n)) $ \(Z :. i :. j) ->
    let a = div i 2
        b = div j 2
        c = mod i 2
        d = mod j 2
     in if abs (a - b) == 1
          then if a == 0
                 then if d == 1
                        then if c == 0
                               then p
                               else 1
                        else 0
                 else if a == n - 1
                        then if d == 0
                               then if c == 0
                                      then 1
                                      else p
                               else 0
                        else if b < a
                               then if d == 0
                                      then if c == 0
                                             then 1
                                             else p
                                      else 0
                               else if d == 0
                                      then 0
                                      else if c == 0
                                             then p
                                             else 1
          else 0

{-# INLINE powerMethodIO #-}
powerMethodIO :: Matrix Double -> Vector Double -> Int -> IO (Vector Double)
powerMethodIO mat vec n = do
  let x = mat #> vec
      maxV = maxElement x
      eigenVec = x / scalar maxV
  -- printf "%d: %f\n" n (rayleighQuotient mat eigenVec)
  -- printf "%03d: %s" n (printList . sumVec . NL.toList $ eigenVec)
  return eigenVec


{-# INLINE powerMethod #-}
powerMethod :: Matrix Double -> Vector Double -> (Double, Vector Double)
powerMethod mat vec =
  let x = mat #> vec
      maxV = maxElement x
      eigenVec = x / scalar maxV
      rq = rayleighQuotient mat eigenVec
   in (rq, eigenVec)

{-# INLINE sumVec #-}
sumVec :: [Double] -> [Double]
sumVec []       = []
sumVec (x:y:xs) = (x + y) : sumVec xs

{-# INLINE printList #-}
printList :: [Double] -> String
printList []     = "\n"
printList (x:xs) = printf "%.3f " x L.++ printList xs

{-# INLINE rayleighQuotient #-}
rayleighQuotient :: Matrix Double -> Vector Double -> Double
rayleighQuotient mat vec = vec <.> (mat #> vec) / vec <.> vec

main = do
  args <- getArgs
  let (nStr:aStr:mStr:_) = args
      n = read nStr :: Int
      a = read aStr :: Double
      m = read mStr :: Int
      mat = tr . ((2 * n) >< (2 * n)) . R.toList . makeMatrix n $ a
      vec = vector . L.take (2 * n) $ [1,1 ..]
      (eigVal', eigVec') = eig mat
      (maxEigVal, maxEigVec) =
        (\(x, y) -> (realPart x, sumVec $ L.map realPart y)) .
        L.maximumBy (\x y -> compare (realPart $ fst x) (realPart $ fst y)) .
        L.zip (NL.toList eigVal') . L.map NL.toList . NL.toColumns $
        eigVec'
      (maxEigVal', maxEigVec') =
        (\(x, y) -> (realPart x, L.map realPart y)) .
        L.maximumBy (\x y -> compare (realPart $ fst x) (realPart $ fst y)) .
        L.zip (NL.toList eigVal') . L.map NL.toList . NL.toColumns $
        eigVec'
      maxV = L.maximumBy (\x y -> compare (abs x) (abs y)) maxEigVec
      (rqs, eigVecs) =
        L.unzip $ L.scanl' (\(_, v) _ -> powerMethod mat v) (0, vec) [1 .. m]
      -- eigVec = L.last eigVecs
      folderPath = "output/test/EndPointMarkovChain"
  createDirectoryIfMissing True folderPath
  eigVec <- M.foldM (powerMethodIO mat) vec [1 .. m]
  -- print . rayleighQuotient mat $ eigVec
  -- putStr . printList . sumVec . NL.toList $ eigVec
  -- print maxEigVal
  -- putStr . printList . L.map (/ maxV) $ maxEigVec
  printf
    "%.3f: %s"
    (L.sum . L.map (/ maxV) $ maxEigVec)
    (printList . L.map (/ maxV) $ maxEigVec)
  plotPathStyle
    [ PNG (folderPath </> "RayleighQuotient.png")
    , Title ("Rayleigh Quotient")
    , XLabel "Iteration"
    ]
    (defaultStyle {plotType = Lines, lineSpec = CustomStyle [LineTitle ""]}) .
    L.zip [1 ..] . L.tail $
    rqs
  let m = 10
      croppedVec =
        L.take (2 * m) . L.drop (div (2 * n - 2 * m) 2) $
        L.map (/ maxV) $ maxEigVec'
  print . L.sum $ croppedVec
