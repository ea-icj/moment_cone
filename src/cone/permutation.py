from .typing import *
from .utils import count, group_by_block

from functools import cached_property
import itertools

class Permutation(tuple[int, ...]): # Remark: hash of p is hash of underlying tuple
    """
    Permutation of S_n represented using the one-line notation.

    The coefficients are the image of range(n) by the permutation.
    So that length computation is faster.
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
        """ Check is permutation is decreasing along each block of given sizes """
        indexes = itertools.accumulate(symmetries, initial=0)
        return all(
            all(a > b for a, b in itertools.pairwise(self[i:j]))
            for i, j in itertools.pairwise(indexes)
        )
    
    def orbit_symmetries(self, symmetries: Iterable[int]) -> Iterable["Permutation"]:
        """
        Permutation inside each block of given sizes
        
        If this is too slow, we may consider the remarks/propositions from:
        - https://stackoverflow.com/questions/19676109/how-to-generate-all-the-permutations-of-a-multiset/
        - https://stackoverflow.com/questions/70057504/speed-up-multiset-permutations
        """
        from sympy.utilities.iterables import multiset_permutations
        indexes = itertools.accumulate(symmetries, initial=0)
        blocks = (multiset_permutations(self[i:j]) for i, j in itertools.pairwise(indexes))
        for p in itertools.product(*blocks):
            yield Permutation(itertools.chain.from_iterable(p))


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
    

class AllPermutationsByLength:
    """
    Catalogue of all permutations of S_n sorted by increasing length

    Since Permutation use a cache for its properties, that means
    this catalogue can use lot of memory to store each permutation,
    it's list of inversion, it's length, ...
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
        
