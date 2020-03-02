module Utils.Parallel
  ( module Control.Parallel.Strategies
  , module Utils.Parallel
  ) where

import           Control.Monad
import           Control.Monad.IO.Class
import           Control.Monad.Parallel       as MP
import           Control.Monad.Trans.Resource
import           Control.Parallel.Strategies
import           Data.Conduit                 as C
import           Data.Conduit.List            as C
import           Data.List                    as L
import           Data.Vector                  as V
import           Prelude                      as P

data ParallelParams = ParallelParams
  { numThread :: Int
  , batchSize :: Int
  } deriving (Show)

{-Parallel Functions-}
{-# INLINE parMapChunk #-}
parMapChunk
  :: ParallelParams -> Strategy b -> (a -> b) -> [a] -> [b]
parMapChunk ParallelParams {numThread = nt} strat f xs =
  (withStrategy (parListChunk (div (P.length xs) nt) strat) . P.map f) xs

{-# INLINE parZipWithChunk #-}
parZipWithChunk :: ParallelParams -> Strategy c -> (a -> b -> c) -> [a] -> [b] -> [c]
parZipWithChunk ParallelParams {numThread = nt} strat f xs =
  withStrategy (parListChunk (div (P.length xs) nt) strat) . P.zipWith f xs

{-# INLINE parZipWith #-}
parZipWith :: Strategy c -> (a -> b -> c) -> [a] -> [b] -> [c]
parZipWith strat f xs  = withStrategy (parList strat) . P.zipWith f xs

{-# INLINE parZipWith3 #-}
parZipWith3
  :: Strategy d -> (a -> b -> c -> d) -> [a] -> [b] -> [c] -> [d]
parZipWith3 strat f xs ys = withStrategy (parList strat) . P.zipWith3 f xs ys

{-# INLINE parZipWith3Chunk #-}
parZipWith3Chunk :: ParallelParams
                 -> Strategy d
                 -> (a -> b -> c -> d)
                 -> [a]
                 -> [b]
                 -> [c]
                 -> [d]
parZipWith3Chunk ParallelParams {numThread = nt} strat f xs ys = withStrategy (parListChunk (div (P.length xs) nt) strat) . P.zipWith3 f xs ys

{- Boxed Vector -}
parMapVector
  :: Strategy b -> (a -> b) -> V.Vector a -> V.Vector b
parMapVector strat f = withStrategy (parTraversable strat) . V.map f

parMapChunkVector :: ParallelParams
                  -> Strategy b
                  -> (a -> b)
                  -> V.Vector a
                  -> V.Vector b
parMapChunkVector ParallelParams{numThread = nt} strat f xs =
  ((withStrategy
      (parVectorChunk (div (V.length xs) nt)
                      strat)) .
   (V.map f)) xs

parZipWithVector :: Strategy c
                 -> (a -> b -> c)
                 -> V.Vector a
                 -> V.Vector b
                 -> V.Vector c
parZipWithVector strat f xs =
  withStrategy (parTraversable strat) . V.zipWith f xs

parZipWithChunkVector :: ParallelParams
                      -> Strategy c
                      -> (a -> b -> c)
                      -> V.Vector a
                      -> V.Vector b
                      -> V.Vector c
parZipWithChunkVector ParallelParams{numThread = nt} strat f xs ys =
  (withStrategy
     (parVectorChunk (div (V.length xs) nt)
                     strat)) .
  (V.zipWith f xs) $
  ys


parZipWith3Vector :: Strategy d
                  -> (a -> b -> c -> d)
                  -> V.Vector a
                  -> V.Vector b
                  -> V.Vector c
                  -> V.Vector d
parZipWith3Vector start f xs ys =
  withStrategy (parTraversable start) . V.zipWith3 f xs ys

parZipWith3ChunkVector :: ParallelParams
                       -> Strategy d
                       -> (a -> b -> c -> d)
                       -> V.Vector a
                       -> V.Vector b
                       -> V.Vector c
                       -> V.Vector d
parZipWith3ChunkVector ParallelParams{numThread = nt} strat f xs ys zs =
  (withStrategy
     (parVectorChunk (div (V.length xs) nt)
                     strat)) .
  (V.zipWith3 f xs ys) $
  zs

parZipWith4ChunkVector :: ParallelParams
                       -> Strategy e
                       -> (a -> b -> c -> d -> e)
                       -> V.Vector a
                       -> V.Vector b
                       -> V.Vector c
                       -> V.Vector d -> V.Vector e
parZipWith4ChunkVector ParallelParams{numThread = nt} strat f xs ys zs as=
  (withStrategy
     (parVectorChunk (div (V.length xs) nt)
                     strat)) .
  (V.zipWith4 f xs ys zs) $
  as

vectorChunk
  :: Int -> V.Vector a -> [V.Vector a]
vectorChunk n vec
  | V.null vec = []
  | otherwise = as : (vectorChunk n bs)
  where
    (as, bs) = V.splitAt n vec

parVectorChunk :: Int -> Strategy a -> Strategy (V.Vector a)
parVectorChunk n strat xs
  | n <= 1 = parTraversable strat xs
  | otherwise =
    fmap V.concat $ parTraversable (evalTraversable strat) (vectorChunk n xs)

parConduit ::
     (NFData b) => ParallelParams -> (a -> b) -> ConduitT a b (ResourceT IO) ()
parConduit parallelParams func = do
  xs <- C.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do sourceList $ parMap rdeepseq func xs
        parConduit parallelParams func)

{-# INLINE parConduitIO #-}
parConduitIO :: ParallelParams -> (a -> IO b) -> ConduitT a b (ResourceT IO) ()
parConduitIO parallelParams func = do
  xs <- C.take (batchSize parallelParams)
  unless
    (L.null xs)
    (do ys <- liftIO $ MP.mapM func xs
        sourceList ys
        parConduitIO parallelParams func)
