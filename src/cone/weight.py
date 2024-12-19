from .typing import *
from .dimension import Dimension

import itertools
import operator

__all__ = (
    "Weight",
)

class Weight:
    """ A weight on which tau can be applied """
    __slots__ = '_weights', 'index'
    _weights: tuple[int, ...]
    index: Optional[int] # Index in the lexicographical order

    def __init__(self, weights: Iterable[int], index: Optional[int] = None):
        self._weights = tuple(weights)
        self.index = index

    @staticmethod
    def all(d: Dimension) -> Iterable["Weight"]:
        """ Returns all possible weights for a given sequence of dimensions, in the lexicographical order """
        for idx, w in enumerate(itertools.product(*(range(di) for di in d))):
            yield Weight(w, idx)

    @staticmethod
    def all_mod_sym_dim(d: Dimension) -> Iterable["Weight"]:
        """
        Returns all decreasing weights modulo the symmetries of d
        
        It assumes that d is arranged so that blocks are contiguous.
        """
        # Using Partition to generate all decreasing weights within a block.
        # The whole weights will be defined as the Cartesian product of the weights for each block.
        from .partition import Partition
        from .utils import group_by_block

        def pad(p: Partition, l: int) -> tuple[int, ...]:
            return p.pad(l)
        
        block_weights = tuple(
            tuple( # Converting to tuple (instead of a map) seems necessary to keep the right hi (FIXME)
                p.pad(hi) # Adding trailing zeros if necessary
                for p in Partition.all_of_height(hi, di - 1)
            )
            for di, hi in group_by_block(d) # Compress returns (value, multiplicity) for each block of d
        )

        for w in itertools.product(*block_weights):
            yield Weight(sum(w, start=())) # Summing tuples is concatenating them

    @staticmethod
    def from_index(d: Dimension, index: int) -> "Weight":
        """ Generate the weight of given index in the lexicographical order """
        weights = []
        tmp_index = index
        stride = tuple(itertools.accumulate(reversed(d[1:]), operator.mul, initial=1))
        for si in reversed(stride):
            weights.append(tmp_index // si)
            tmp_index -= weights[-1] * si
        return Weight(weights, index)

    def index_in(self, d: Dimension, use_internal_index: bool = True) -> int:
        """
        Returns index of this weight in the lexicographical order for given dimensions (see `all` method)
        
        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """
        if not use_internal_index or self.index is None:
            stride = itertools.accumulate(reversed(d[1:]), operator.mul, initial=1)
            self.index = sum(v * s for v, s in zip(reversed(self._weights), stride))
        return self.index

    def orbit_symmetries(self, symmetries: Iterable[int]) -> Iterable["Weight"]:
        """
        Permutation inside each block of given sizes

        Example:
        >>> p = Weight((2, 2, 4, 1, 2, 1, 4))
        >>> orbits = p.orbit_symmetries((3, 3, 1))
        >>> for chi in orbits:
        ...     print(chi)
        Weight((2, 2, 4, 1, 1, 2, 4))
        Weight((2, 2, 4, 1, 2, 1, 4))
        Weight((2, 2, 4, 2, 1, 1, 4))
        Weight((2, 4, 2, 1, 1, 2, 4))
        Weight((2, 4, 2, 1, 2, 1, 4))
        Weight((2, 4, 2, 2, 1, 1, 4))
        Weight((4, 2, 2, 1, 1, 2, 4))
        Weight((4, 2, 2, 1, 2, 1, 4))
        Weight((4, 2, 2, 2, 1, 1, 4))
        """
        from .utils import orbit_symmetries
        for chi in orbit_symmetries(self._weights, symmetries):
            yield Weight(chi)

    def __len__(self) -> int:
        return len(self._weights)

    @overload
    def __getitem__(self, idx: int) -> int:
        ...

    @overload
    def __getitem__(self, idx: slice) -> tuple[int, ...]:
        ...

    def __getitem__(self, idx: int | slice) -> int | tuple[int, ...]:
        return self._weights[idx]
    
    def __iter__(self) -> Iterator[int]:
        return iter(self._weights)

    def __le__(self, other: "Weight") -> bool:
        """ Implementation of self <= other (partial ordering)"""
        return all(ws >= wo for ws, wo in zip(self, other))

    def __ge__(self, other: "Weight") -> bool:
        """ Implementation of self >= other (partial ordering)"""
        return all(ws <= wo for ws, wo in zip(self, other))
    
    def __eq__(self, other: object) -> bool:
        """ Equality between two weights (ignoring index) """
        if not isinstance(other, Weight):
            return NotImplemented
        return self._weights == other._weights
    
    def __repr__(self) -> str:
        return f"Weight({self._weights}" + (f", idx: {self.index}" if self.index is not None else "") + ")"

