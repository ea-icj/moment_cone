from .typing import *
from functools import cached_property
import itertools

class Permutation(tuple[int, ...]):
    """
    Permutation represented using the one-line notation.

    So that length computation is faster.
    """

    @cached_property
    def length(self) -> int:
        """ Length of the permutation, ie number of inversions """
        return sum(
            1 if self[i] > self[j] else 0
            for i, j in itertools.combinations(range(len(self)), 2)
        )
    
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
    