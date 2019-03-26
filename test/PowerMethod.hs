module PowerMethod where

import           Control.Monad         as M
import           Data.Complex
import           Data.List             as L
import           Numeric.LinearAlgebra as NL
import           System.Environment
import           System.Random
import           Text.Printf

{-# INLINE iterateM #-}
iterateM :: Int -> (b -> a -> IO a) -> b -> a -> IO a
iterateM 0 _ _ x = return x
iterateM n f x y = do
  z <- f x y
  iterateM (n - 1) f x $! z
  

  
{-# INLINE getComplexStr #-}
getComplexStr :: (Complex Double) -> String
getComplexStr (a :+ b) = printf "%.2f :+ %.2f" a b

{-# INLINE getComplexStr' #-}
getComplexStr' :: (Complex Double) -> String
getComplexStr' x =
  let (a, b) = polar x
   in printf "(%.2f, %.2f)" a b

{-# INLINE getListStr #-}
getListStr :: (a -> String) -> [a] -> String
getListStr f (x:xs) =
  printf "[%s%s]" (f x) . L.concatMap (\x -> " ," L.++ f x) $ xs

main = do
  (nStr:_) <- getArgs
  let range = (-1, 1)
  matInit1 <- M.replicateM 4 (randomRIO range) :: IO [Double]
  matInit2 <- M.replicateM 4 (randomRIO range) :: IO [Double]
  vecInit1 <- M.replicateM 2 (randomRIO range) :: IO [Double]
  vecInit2 <- M.replicateM 2 (randomRIO range) :: IO [Double]
  let n = read nStr :: Int
      mat' = (2 >< 2) $ L.zipWith (:+) matInit1 matInit2
      vec' = NL.fromList $ L.zipWith (:+) vecInit1 vecInit2
      (eigVal, eigVec) = eig mat'
      eigVecs = toColumns eigVec
  print mat'
  out <-
    iterateM
      n
      (\mat vec -> do
         let newVec = mat #> vec
             maxMag = L.maximum . L.map magnitude . NL.toList $ newVec
             (a:b:_) = NL.toList newVec
             phaseNorm x =
               let (m, p) = polar x
                in mkPolar 1 p
             norm =
               if magnitude a > magnitude b
                 then a
                 else b
             normalizedNewVec = newVec / (scalar norm)
         printf
           "Maximum magnitude: %.2f %s\n"
           maxMag
           (getListStr getComplexStr' . NL.toList $ normalizedNewVec)
         return normalizedNewVec)
      mat'
      vec'
  let s = sqrt . L.sum . L.map (\x -> (magnitude x) ^ 2) . NL.toList $ out
  printf
    "%s\n"
    (getListStr getComplexStr' . L.map (/ (s :+ 0)) . NL.toList $ out)
  M.zipWithM_
    (\val vec ->
       printf
         "%.2f %s\n"
         (magnitude val)
         (getListStr getComplexStr' . NL.toList $ vec))
    (NL.toList eigVal)
    eigVecs
