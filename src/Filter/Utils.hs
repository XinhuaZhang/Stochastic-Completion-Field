{-# LANGUAGE TypeOperators #-}
module Filter.Utils where

import           Data.Array.Repa as R
import           Data.List       as L

{-# INLINE makeFilter #-}
makeFilter ::
     (Source r e, Shape sh)
  => R.Array r (sh :. Int :. Int) e
  -> R.Array D (sh :. Int :. Int) e
makeFilter arr =
  let (cols:rows:_) = L.take 2 . listOfShape . extent $ arr
   in R.backpermute
        (extent arr)
        (\(sh :. i :. j) ->
           let halfRows = div rows 2
               halfCols = div cols 2
               x =
                 if i < halfRows
                   then i + halfRows
                   else i - halfRows
               y =
                 if j < halfCols
                   then j + halfCols
                   else j - halfCols
            in (sh :. x :. y))
        arr
