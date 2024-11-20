"""
Comparison of different implementation to generate all partitions
of S_n and compute their length, inversions and the inverse.

There are three different versions compared here:
1. the original code that relies on Sage,
2. the new code with permutations represented as one-line (image of range(n)),
3. a variation of this new code where each permutation is associated to one unique instance
   so that to avoid storing and computing permutations (and there attributes) multiple times.
   Eg when computing the inverse, the inverse of the inverse is automatically linked to the
   current permutation so that when instantiating this inverse elsewhere, it's inverse is
   already available.

In addition, the original comparison code was also computing the inversions of the inverse
due to a misunderstanding of what the Sage code do (it compute the inversions of the inverse
so that to get the inversions in the index format). In that case, the 3rd code was two times
faster than the 2nd code. Without this extra computation, code duration of this two variations
are of the same order.

However, the third version can reduce memory consumption.

Here are some computation time results (in ms) :

| S_n |   sage   | nocache | nocache2 |   cache   |
| --- | -------- | ------- | -------- | --------- |
|  4  |    0.921 |   0.077 |   0.086  |   0.086ms |
|  5  |    4.502 |   0.427 |   0.450  |   0.455ms |
|  6  |   34.676 |   2.946 |   3.117  |   3.089ms |
|  7  |  640.326 |  25.519 |  25.663  |  25.891ms |
|  8  | 5332.208 | 234.243 | 228.455  | 224.559ms |

"""

from pathlib import Path
import sys
SOURCE_PATH = Path.cwd() / 'src'
sys.path.append(str(SOURCE_PATH))

from sage.all import SymmetricGroup, Permutation as SagePermutation

from cone.typing import *
from cone.utils import count, group_by_block

from functools import cached_property
import itertools

def method_sage(n: int):
    """ Initial way using Sage """
    Sd = SymmetricGroup(n)
    data = []
    for w in Sd:
        invers = [list(beta.cycle_tuples()[0]) for beta in  w.inverse().inversions()]
        data.append([
            w.length(),
            SagePermutation(w),
            invers,
            SagePermutation(w.inverse())
        ])
    return data


class Permutation(tuple[int, ...]): # Remark: hash of p is hash of underlying tuple
    """
    Permutation represented using the one-line notation.

    So that length computation is faster.
    """
    @cached_property
    def inversions(self) -> tuple[tuple[int, int], ...]:
        """ Sequence of the indexes of all the inversions """
        return tuple(filter(
            lambda ij: self[ij[0]] > self[ij[1]],
            itertools.combinations(range(len(self)), 2)
        ))
    
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

    @staticmethod
    def all(n: int) -> Iterable["Permutation"]:
        """ Returns all the possible permutation of S_n """
        for p in itertools.permutations(range(n)):
            yield Permutation(p)

def method_nocache(n: int):
    """ Using one-line notation of the permutations """
    data = tuple(Permutation.all(n))
    for p in data:
        p.length, p.inverse, p.inversions
    return data


class PermutationCache(tuple[int, ...]): # Remark: hash of p is hash of underlying tuple
    """
    Permutation represented using the one-line notation.

    So that length computation is faster.
    """
    all_permutations: dict["PermutationCache", "PermutationCache"] = {}
    
    def __new__(cls, coefficients):
        p = super().__new__(cls, coefficients)
        return cls.all_permutations.setdefault(p, p)

    @cached_property
    def inversions(self) -> tuple[tuple[int, int], ...]:
        """ Sequence of the indexes of all the inversions """
        return tuple(filter(
            lambda ij: self[ij[0]] > self[ij[1]],
            itertools.combinations(range(len(self)), 2)
        ))
    
    @cached_property
    def length(self) -> int:
        """ Length of the permutation, ie number of inversions """
        return count(self.inversions)
    
    @cached_property
    def inverse(self) -> "PermutationCache":
        """ Inverse of the permutation """
        inv = [0] * len(self)
        for i, pi in enumerate(self):
            inv[pi] = i
        p_inv = PermutationCache(inv)
        p_inv.inverse = self
        return p_inv

    @staticmethod
    def all(n: int) -> Iterable["PermutationCache"]:
        """ Returns all the possible permutation of S_n """
        for p in itertools.permutations(range(n)):
            yield PermutationCache(p)

def method_cache(n: int):
    """ Using a cache of the one-line notation of the permutations """
    PermutationCache.all_permutations.clear()
    data = tuple(PermutationCache.all(n))
    for p in data:
        p.length, p.inverse, p.inversions
    return data

class Permutation2(tuple[int, ...]): # Remark: hash of p is hash of underlying tuple
    """
    Permutation represented using the one-line notation.

    Small optimization where the inverse of the inverse is already filled.
    (it will not speed up this comparison)
    """
    @cached_property
    def inversions(self) -> tuple[tuple[int, int], ...]:
        """ Sequence of the indexes of all the inversions """
        return tuple(filter(
            lambda ij: self[ij[0]] > self[ij[1]],
            itertools.combinations(range(len(self)), 2)
        ))
    
    @cached_property
    def length(self) -> int:
        """ Length of the permutation, ie number of inversions """
        return count(self.inversions)
    
    @cached_property
    def inverse(self) -> "Permutation2":
        """ Inverse of the permutation """
        inv = [0] * len(self)
        for i, pi in enumerate(self):
            inv[pi] = i
        p_inv = Permutation2(inv)
        p_inv.inverse = self
        return p_inv

    @staticmethod
    def all(n: int) -> Iterable["Permutation2"]:
        """ Returns all the possible permutation of S_n """
        for p in itertools.permutations(range(n)):
            yield Permutation2(p)

def method_nocache2(n: int):
    """ Using one-line notation of the permutations """
    data = tuple(Permutation2.all(n))
    for p in data:
        p.length, p.inverse, p.inversions
    return data


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        "Comparing methods to compute informations for all permutations of S_n",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("n", type=int, help="The dimension of S_n")
    parser.add_argument("--methods", choices=["sage", "nocache", "nocache2" ,"cache"], nargs="*", default=["sage", "nocache", "nocache2", "cache"], help="Methods to test")
    config = parser.parse_args()

    print(f"Computing permutation in S_{config.n}:")

    from timeit import Timer
    for method in config.methods:
        print(f"\tusing method {method:<7}: ", end='', flush=True)
        timer = Timer(
            stmt=f"method_{method}({config.n})",
            setup=f"from {__name__} import method_{method}"
        )
        cnt, duration = timer.autorange()
        print(f"{1000 * duration / cnt:9.3f}ms ({cnt} runs)")


