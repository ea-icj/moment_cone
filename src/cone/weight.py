from .typing import *
import itertools
import operator

"""
TODO:
- OrbitWeightSd method (needs multiset_permutations)
"""

class Weight(Sequence[int]):
    """ A weight on which tau can be applied """
    __slots__ = '_weights', 'index'
    _weights: tuple[int, ...]
    index: Optional[int] # Index in the lexicographical order

    def __init__(self, weights: Iterable[int], index: Optional[int] = None):
        self._weights = tuple(weights)
        self.index = index

    @staticmethod
    def all(d: Iterable[int]) -> Iterable["Weight"]:
        """ Returns all possible weights for a given sequence of dimensions, in the lexicographical order """
        for idx, w in enumerate(itertools.product(*(range(di) for di in d))):
            yield Weight(w, idx)

    @staticmethod
    def from_index(d: Sequence[int], index: int) -> "Weight":
        """ Generate the weight of given index in the lexicographical order """
        weights = []
        tmp_index = index
        stride = tuple(itertools.accumulate(reversed(d[1:]), operator.mul, initial=1))
        for si in reversed(stride):
            weights.append(tmp_index // si)
            tmp_index -= weights[-1] * si
        return Weight(weights, index)

    def index_in(self, d: Sequence[int]) -> int:
        """ Returns index of this weight in the lexicographical order for given dimensions (see `all` method)"""
        stride = itertools.accumulate(reversed(d[1:]), operator.mul, initial=1)
        return sum(v * s for v, s in zip(reversed(self._weights), stride))
    
    def __len__(self) -> int:
        return len(self._weights)
    
    def __getitem__(self, idx: int) -> int:
        return self._weights[idx]
    
    def __iter__(self) -> Iterator[int]:
        return iter(self._weights)

    def __le__(self, other: "Weight") -> bool:
        """ Implementation of self <= other (partial ordering)"""
        return all(ws >= wo for ws, wo in zip(self, other))

    def __ge__(self, other: "Weight") -> bool:
        """ Implementation of self >= other (partial ordering)"""
        return all(ws <= wo for ws, wo in zip(self, other))
    
    def __eq__(self, other: "Weight") -> bool:
        """ Equality between two weights (ignoring index) """
        return self._weights == other._weights

