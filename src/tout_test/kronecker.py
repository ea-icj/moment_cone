__all__ = (
    'KroneckerWeight',
    'KroneckerRepresentation',
)

from functools import cached_property
import itertools

from .typing import *
from .representation import Representation
from .weight import Weight as WeightBase
from .linear_group import LinearGroup
from .rings import Vector, vector

class KroneckerWeight(WeightBase):
    """
    Weight class for the Kronecker representation
    """
    __weights: tuple[int, ...]

    def __init__(self, G: LinearGroup, weights: Iterable[int], **kwargs):
        super().__init__(G, **kwargs)
        self.__weights = tuple(weights)
        assert (
            len(self.__weights) == len(self.G)
            and all(0 <= w < g for w, g in zip(self.__weights, self.G))
        ), "Invalid weight"

    @cached_property
    def as_vector(self) -> Vector:
        from .rings import ZZ
        v = vector(ZZ, self.G.rank)
        for shift, x in zip(itertools.accumulate(self.G, initial=0), self.__weights):
            v[shift + x] = 1
        return v

    def __le__(self, other: WeightBase) -> bool:
        """
        Implementation of self <= other (partial ordering)
        
        Example:
        >>> G = LinearGroup((3, 2))
        >>> K = KroneckerRepresentation(G)
        >>> chi1 = K.weight((2, 1))
        >>> chi2 = K.weight((2, 0))
        >>> chi1 <= chi2
        True

        >>> from .weight import Weight
        >>> Weight(G, as_vector=chi1.as_vector) <= Weight(G, as_vector=chi2.as_vector)
        True
        """
        if not isinstance(other, KroneckerWeight):
            return NotImplemented
        return all(ws >= wo for ws, wo in zip(self, other))
    
    def leq(self,
            other: "WeightBase",
            symmetries: Optional[Iterable[int]] = None) -> bool:
        if symmetries is None:
            return self <= other
        else:
            return super().leq(other, symmetries)
        
    def __eq__ (self, other: object) -> bool:
        if not isinstance(other, KroneckerWeight):
            return NotImplemented
        return self.G == other.G and self.__weights == other.__weights

    def __repr__(self) -> str:
        return f"KroneckerWeight({self.__weights}" + (f", idx: {self.index}" if self.index is not None else "") + ")"

    def __len__(self) -> int:
        return len(self.__weights)

    @overload
    def __getitem__(self, idx: int) -> int:
        ...

    @overload
    def __getitem__(self, idx: slice) -> tuple[int, ...]:
        ...

    def __getitem__(self, idx: int | slice) -> int | tuple[int, ...]:
        return self.__weights[idx]
    
    def __iter__(self) -> Iterator[int]:
        return iter(self.__weights)

class KroneckerRepresentation(Representation):
    Weight = KroneckerWeight
    
    @cached_property
    def dim_cone(self) -> int:
        return self.G.rank - len(self.G) + 1
    
    @cached_property
    def dim(self) -> int:
        from .utils import prod
        return prod(self.G)
    
    def all_weights(self) -> Iterable[WeightBase]:
        """
        Returns all possible weights for a given sequence of dimensions, in the lexicographical order
        """
        for idx, w in enumerate(itertools.product(*(range(di) for di in self.G))):
            yield KroneckerWeight(self.G, w, index=idx)

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        """
        Returns index of the weight in the lexicographical order for the current representation (see `all_weights` method).
        
        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """
        if not use_internal_index or chi.index is None:
            if not isinstance(chi, KroneckerWeight):
                raise ValueError("Invalid weight representation")        
            from operator import mul
            stride = itertools.accumulate(reversed(self.G[1:]), mul, initial=1)
            chi.index = sum(v * s for v, s in zip(reversed(chi), stride))
        return chi.index
