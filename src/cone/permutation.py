from .typing import *
from .utils import count, group_by_block
from .blocks import Blocks

from functools import cached_property
import itertools

__all__ = (
    "Permutation",
    "AllPermutationsByLength",
)

class Permutation(tuple[int, ...]): # Remark: hash of p is hash of underlying tuple
    """
    Permutation of S_n represented using the one-line notation.

    The coefficients are the image of range(n) by the permutation.
    So that length computation is faster.

    Example:
    >>> from cone.permutation import Permutation
    >>> p = Permutation((3, 2, 1, 4, 0, 5))
    >>> p
    Permutation((3, 2, 1, 4, 0, 5))
    >>> p.n
    6
    >>> p.inverse
    Permutation((4, 2, 1, 0, 3, 5))
    >>> p(p.inverse)
    (0, 1, 2, 3, 4, 5)
    >>> p.length
    7
    >>> for p in Permutation.all(3):
    ...     print(p)
    Permutation((0, 1, 2))
    Permutation((0, 2, 1))
    Permutation((1, 0, 2))
    Permutation((1, 2, 0))
    Permutation((2, 0, 1))
    Permutation((2, 1, 0))
    """
    @property
    def n(self) -> int:
        return len(self)
    
    @cached_property
    def inversions(self) -> tuple[tuple[int, int], ...]:
        """ Sequence of the indexes of all the inversions """
        return tuple(filter(
            lambda ij: self[ij[0]] > self[ij[1]],
            itertools.combinations(range(self.n), 2)
        ))
    
    @cached_property
    def length(self) -> int:
        """ Length of the permutation, ie number of inversions """
        return count(self.inversions)
    
    @cached_property
    def inverse(self) -> "Permutation":
        """ Inverse of the permutation """
        inv = [0] * self.n
        for i, pi in enumerate(self):
            inv[pi] = i
        p_inv = Permutation(inv)
        p_inv.inverse = self
        return p_inv

    def __call__(self, s: Sequence[T]) -> tuple[T, ...]:
        """ Apply the permutation to a given sequence """
        return tuple(s[pi] for pi in self)
    
    def __repr__(self) -> str:
        return f"Permutation({super().__repr__()})"
    
    def is_min_rep(self, symmetries: Iterable[int]) -> bool:
        """ Check if permutation is decreasing along each block of given sizes """
        return all(
            all(a > b for a, b in itertools.pairwise(block))
            for block in Blocks.from_flatten(self, symmetries)
        )

    @staticmethod
    def all(n: int) -> Iterable["Permutation"]:
        """ Returns all the possible permutation of S_n """
        for p in itertools.permutations(range(n)):
            yield Permutation(p)

    @staticmethod
    def all_of_length(n: int, l: int) -> Iterable["Permutation"]:
        """
        Returns all permutations of S_n with given length l

        Better use AllPermutationsByLength for a repeated call with difference lengths
        """
        # More efficient way ?
        return filter(lambda p: p.length == l, Permutation.all(n))
    
    @staticmethod
    def from_cycles(n: int, *cycles: Sequence[int]) -> "Permutation":
        """
        Returns a permutation from a list of cycles
        
        >>> Permutation.from_cycles(5, (2, 1), (3, 0, 4))
        Permutation((4, 2, 1, 0, 3))
        """
        p = list(range(n))
        for cycle in cycles:
            for i, j in itertools.pairwise(cycle):
                p[i] = j
            p[cycle[-1]] = cycle[0]
        return Permutation(p)
    
    def transpose(self, i: int, j: int) -> "Permutation":
        """
        Apply a transposition between ith and jth positions

        >>> p = Permutation((4, 2, 1, 0, 3))
        >>> p.transpose(3, 1)
        Permutation((4, 0, 1, 2, 3))
        """
        if i > j:
            i, j = j, i
        return Permutation(
            self[:i] + self[j:j+1] + self[i+1:j] + self[i:i+1] + self[j+1:]
        )

    @cached_property
    def covering_relations_strong_Bruhat(self) -> tuple["Permutation", ...]:
        """
        Covering relations strong Bruhat

        Liste of v <= self for the ??? order so that v <= self is wrong
        for the weak order and v has a length that is equal to the length of self minus 1.
        """
        n = len(self)
        return tuple(
            p
            for i, j in itertools.combinations(range(n), 2)
            for p in (self.transpose(i, j),)
            if p.length == self.length - 1
        )


    @staticmethod
    def all_transpositions(n: int) -> Iterable["Permutation"]:
        """
        Returns all transpositions
        
        >>> for p in Permutation.all_transpositions(3):
        ...     print(p)
        Permutation((1, 0, 2))
        Permutation((2, 1, 0))
        Permutation((0, 2, 1))
        """
        return (Permutation.from_cycles(n, c) for c in itertools.combinations(range(n), 2))


class AllPermutationsByLength:
    """
    Catalogue of all permutations of S_n sorted by increasing length

    Since Permutation use a cache for its properties, that means
    this catalogue can use lot of memory to store each permutation,
    it's list of inversion, it's length, ...

    Examples:
    >>> ap = AllPermutationsByLength(3)
    >>> len(ap) # Number of different lengths
    4
    >>> ap.max_length # Maximal length (should be len(ap) - 1)
    3
    >>> ap[1] # Permutations of length 1
    (Permutation((0, 2, 1)), Permutation((1, 0, 2)))
    >>> ap[:2] # All permutations of length <= 1
    (Permutation((0, 1, 2)), Permutation((0, 2, 1)), Permutation((1, 0, 2)))
    """
    __slots__ = 'permutations', 'indexes'
    permutations: tuple[Permutation, ...]
    indexes: tuple[int, ...]

    def __init__(self, n: int):
        self.permutations = tuple(sorted(
            Permutation.all(n),
            key=lambda p: p.length
        ))

        # Storing index when each length begins
        # Note that it assumes that there is not hole in the lengths (must be the case)
        self.indexes = tuple(itertools.accumulate(
            (size for value, size in group_by_block(p.length for p in self.permutations)),
            initial=0
        ))

    def __len__(self) -> int:
        """
        Returns the number of different lengths (from 0 to maximal length included)
        """
        return len(self.indexes) - 1
    
    @property
    def max_length(self) -> int:
        """ Returns the maximal length of a permutation """
        return len(self) - 1
    
    def __getitem__(self, idx: int | slice) -> tuple[Permutation, ...]:
        """ Access permutations by length or range of length """
        def fix_neg(i: int) -> int:
            return len(self) + i if i < 0 else i
        
        if isinstance(idx, slice):
            assert idx.step in (None, 1), "Only contiguous slice are supported"
            start = 0 if idx.start is None else fix_neg(idx.start)
            stop = len(self) if idx.stop is None else fix_neg(idx.stop)
            return self.permutations[self.indexes[start]:self.indexes[stop]]
        else:
            idx = fix_neg(idx)
            return self.permutations[self.indexes[idx]:self.indexes[idx + 1]]
        
