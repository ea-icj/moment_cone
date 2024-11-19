from .typing import *
from .utils import count

from functools import cached_property
import itertools

class Permutation(tuple[int, ...]):
    """
    Permutation represented using the one-line notation.

    So that length computation is faster.
    """
    @property
    def inversions(self) -> Iterable[tuple[int, int]]:
        """ Sequence of the indexes of all the inversions """
        return filter(
            lambda ij: self[ij[0]] > self[ij[1]],
            itertools.combinations(range(len(self)), 2)
        )
    
    @cached_property
    def length(self) -> int:
        """ Length of the permutation, ie number of inversions """
        return count(self.inversions)
    
    @cached_property
    def inverse(self) -> "Permutation":
        """ Inverse of the permutation """
        inv = [0] * len(self)
        for i, pi in enumerate(self):
            inv[pi] = i
        return Permutation(inv)

    def __call__(self, s: Sequence[T]) -> tuple[T, ...]:
        """ Apply the permutation to a given sequence """
        return tuple(s[pi] for pi in self)
    
    def __repr__(self) -> str:
        return f"Permutation({super().__repr__()})"
    
    @staticmethod
    def all(n: int) -> Iterable["Permutation"]:
        """ Returns all the possible permutation of S_n """
        for p in itertools.permutations(range(n)):
            yield Permutation(p)

    @staticmethod
    def all_of_length(n: int, l: int) -> Iterable["Permutation"]:
        """ Returns all permutations of S_n with given length l """
        # More efficient way ?
        return filter(lambda p: p.length == l, Permutation.all(n))
    