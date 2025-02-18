from .typing import *
from .utils import count, group_by_block, expand_blocks, is_increasing, multiset_permutations
from .blocks import Blocks
from .group import *

from functools import cached_property
import itertools

__all__ = (
    "Permutation",
    "AllPermutationsByLength",
)

class Permutation(tuple[int, ...]): # Remark: hash of p is hash of underlying tuple
    """
    Permutation of S_n represented using the one-line notation.

    The coefficients are the image of range(n) by the permutation, ie `p(range(n)) == p`.
    
    Using that convention, length computation is faster.

    Example:
    >>> from cone.permutation import Permutation
    >>> p = Permutation((3, 2, 1, 4, 0, 5))
    >>> p
    Permutation((3, 2, 1, 4, 0, 5))
    >>> p(range(6))
    (3, 2, 1, 4, 0, 5)
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
        """
        Check if permutation is increasing along each block of given sizes

        Lie theoretically speaking: 
        symmetries encode the size of blocks of a Levi L, 
        Denote its Weyl group by W_L
        w satisfies is_min_rep if it is a representative of W/W_L of minimal length

        Example:
        >>> p = Permutation((1, 2, 3, 3, 4, 2, 3, 4, 5))
        >>> p.is_min_rep((3, 2, 4))
        True
        >>> p.is_min_rep((2, 3, 4))
        False
        """
        return all(
            all(a < b for a, b in itertools.pairwise(block))
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
    def all_min_rep(symmetries: Iterable[int]) -> Iterable["Permutation"]:
        """
        Returns all permutations of S_n that are increasing along each block of given symmetries.

        n is defined by the sum of the block sizes.

        Examples:
        >>> n = 7
        >>> symmetries = (2, 3, 2)
        >>> naive_way = filter(
        ...     lambda p: p.is_min_rep(symmetries),
        ...     Permutation.all(n)
        ... )
        >>> this_way = Permutation.all_min_rep(symmetries)
        >>> all(
        ...     a == b
        ...     for a, b in itertools.zip_longest(naive_way, this_way, fillvalue=None)
        ... )
        True
        """
        def block_recurs(
                seq: tuple[int, ...],
                length: int,
                start: int = 0
            ) -> Generator[tuple[tuple[int, ...], tuple[int, ...]]]:
            """ For one block of symmetry of given length, returns all strictly
            increasing permutations of seq (head) alongside the remaining
            values to permute (tail)"""
            if length == 0:
                yield (), seq[start:]
            else:
                for i in range(start, len(seq) + 1 - length):
                    for head, tail in block_recurs(seq, length - 1, start=i+1):
                        yield seq[i:i+1] + head, seq[start:i] + tail

        def all_recurs(
                seq: tuple[int, ...],
                sym_tail: tuple[int, ...]
            ) -> Generator[tuple[tuple[int, ...], ...]]:
            """ For a given sequence and symmetries, returns all strictly
            increasing permutations of the sequence by block of symmetry
            as a tuple of one permutation per block """
            if len(sym_tail) == 0:
                yield (),
            else:
                for head, rem_seq in block_recurs(seq, sym_tail[0]):
                    for tail in all_recurs(rem_seq, sym_tail[1:]):
                        yield head, *tail

        symmetries = tuple(symmetries)
        n = sum(symmetries)
        for parts in all_recurs(tuple(range(n)), symmetries):
            yield Permutation(itertools.chain.from_iterable(parts))

        
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
    
    @staticmethod
    def embeddings_mod_sym(G: LinGroup, Gred: LinGroup)-> Iterable["Permutation"]:
        """
        List of permutations of e that are at most d

        d and e are list of integers of the same length  (typically, dimensions), each in decreasing order.

        Returns the list of permutation of e (each encoded by a permutation of the indices) such that the value in the i-th component of the permuted e is at most d[i]
        Outputs are irredundant modulo symmetries of e and d
        
        Example:
        >>> d = [4, 4, 3, 3, 2]
        >>> e = [4, 3, 3, 2, 1]
        >>> emb = list(Permutation.embeddings_mod_sym(d, e))
        >>> for pe in emb:
        ...     print(pe)
        Permutation((0, 1, 2, 3, 4))
        Permutation((0, 1, 2, 4, 3))
        Permutation((0, 3, 1, 2, 4))
        Permutation((0, 4, 1, 2, 3))
        >>> for pe in emb:
        ...     print(pe(e))
        (4, 3, 3, 2, 1)
        (4, 3, 3, 1, 2)
        (4, 2, 3, 3, 1)
        (4, 1, 3, 3, 2)
        """
        # TODO: check if group_by_block and expand_blocs can be rewritten in another way
        eg = list(group_by_block(Gred))
        dg = list(group_by_block(G))
        partial_sum_mult_d = [0] + list(itertools.accumulate([x for _,x in dg])) #indexes of the first element in each block
        partial_sum_mult_e = [0] + list(itertools.accumulate([x for _,x in eg]))
        indices_eg = expand_blocks([i for i, x in enumerate(eg)], [x[1] for i, x in enumerate(eg)]) # same shape as Gred, but with 0,1,2 as values exple : [4,3,3,1] --> [0,1,1,2]
        for ep in multiset_permutations(indices_eg): # Any permutation of the list Gred
            p_i = partial_sum_mult_e[:-1] # Same as length of Gred

            indices_e = []
            for i in range(len(ep)):
                indices_e.append(p_i[ep[i]]) # p_i[ep[i]] index in Gred of a value of the block indexed by ep[i]
                p_i[ep[i]] += 1 # we add 1 to take the next value the next time

            if all(Gred[indices_e[i]] <= G[i] for i in range(len(ep))) and all(is_increasing(indices_e[a:b]) for a,b in itertools.pairwise(partial_sum_mult_d)):
                # The first all check that the permutted Gred is still at most G
                # The second all check that in a block of G we took elements of Gred  from left to right                
                yield Permutation(indices_e)
    
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
    def covering_relations_strong_Bruhat_old(self) -> tuple["Permutation", ...]:
        """
        Covering relations strong Bruhat

        Liste of v <= self for the Bruhat order so that v <= self is wrong
        for the weak Bruhat order and v has a length that is equal to the length of self minus 1.
        """
        n = len(self)
        return tuple(
            p
            for i, j in itertools.combinations(range(n), 2)
            for p in (self.transpose(i, j),)
            if p.length == self.length - 1
        )

    @cached_property
    def covering_relations_strong_Bruhat(self) -> tuple["Permutation", ...]:
        """
        Covering relations strong Bruhat

        Liste of v <= self for the Bruhat order so that v <= self is wrong
        for the weak Bruhat order and v has a length that is equal to the length of self minus 1.
        """
        n = len(self)
        return tuple(
            p
            for i, j in itertools.combinations(range(n), 2)
            for p in [Permutation(list(self[:i]) + [self[j]] + list(self[i+1:j]) + [self[i]] + list(self[j+1:]))]
            if (self[i] > self[j]+1) and (p.length == self.length - 1)
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

    It is worth noting that this class ensure uniqueness of an instance for a given n
    so that it can be constructed multiple times for a same n without a computation time penalty.
    >>> ap2 = AllPermutationsByLength(3)
    >>> ap is ap2
    True
    """
    all_instances: dict[int, "AllPermutationsByLength"] = {}

    __slots__ = 'permutations', 'indexes'
    permutations: tuple[Permutation, ...]
    indexes: tuple[int, ...]

    def __new__(cls, n: int):
        """ Construction with reusing of already computed permutations for given n """
        try:
            return cls.all_instances[n]
        except KeyError:
            pass

        self = super().__new__(cls)
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

        cls.all_instances[n] = self
        return self

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
        
