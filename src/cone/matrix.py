from .typing import *

"""
TODO: it seems that this partial matrix representation is not really needed.
It should be more useful to have a more clear representation of a sequence as multiple blocks (column here) depending on a given list of sizes of the blocks.
This thing is used everywhere, not only for tau (look at the usage of group_by_block).

We could eventually keep a mutable version like the one below and a freezed one
that could be optimized as a one tuple like in the original code.

""" 
class PartialMatrix(Generic[T]):
    """ Kind of sparse matrix with fixed number of column, each of them containing a variable number of values """
    __slots__ = '_data', 'shape'
    _data: tuple[list[T], ...] # Internal representation of the columns
    shape: tuple[int, int] # Shape of the array (maximal number of rows and fixed number of cols)

    def __init__(self, nrows: int, ncols: int):
        """ Initialize empty partial matrix of given shape """
        self._data = tuple([] for _ in range(ncols))
        self.shape = nrows, ncols
    
    def size(self, idx: int) -> int:
        """ Number of rows for the given columns """
        return len(self._data[idx])
    
    @property
    def sizes(self) -> tuple[int, ...]:
        """ Number of rows for all columns """
        return tuple(len(col) for col in self._data)
    
    def __getitem__(self, idx: tuple[int, int]) -> T:
        """ Read-only access to an element using classical (i,j) order """
        i, j = idx
        return self._data[j][i]
    
    def column(self, j: int) -> list[T]:
        """ Access to a given column """
        return self._data[j]
    
    def append(self, j: int, value: T) -> None:
        """ Append a value to the given column """
        self._data[j].append(value)

    def extend(self, j: int, values: Iterable[T]) -> None:
        """ Extend a column with given values """
        self._data[j].extend(values)

    @property
    def columns(self) -> tuple[list[T], ...]:
        """ Access to all columns """
        return self._data

    def __repr__(self) -> str:
        return " | ".join(" ".join(map(str, c)) for c in self.columns)


# TODO: un genre de FrozenPartialMatrix non mutable et mieux optimisé (on pourrait reprendre l'idée de tout mettre en ligne)
# En plus, ça permettrait de traiter également les nombreux découpages par blocs de symétrie qui interviennent dans le code.
# D'ailleurs, faudrait se reposer la question de l'intérêt de cette "matrice partielle", utilise-t-on vraiment l'indiçage i,j ? A-t-on besoin de modifier les colonnes ? etc.

