__all__ = (
    'Weight',
    'WeightAsList',
    'WeightAsListOfList',
)

from functools import cached_property
import itertools

from .typing import *
from .linear_group import LinearGroup

if TYPE_CHECKING:
    from .rings import Vector


class Weight:
    """ A weight on which tau can be applied """
    G: LinearGroup
    index: Optional[int] #: Index for an optional generation order
    mult: int #: Multiplicity of this weight

    def __init__(self,
                 G: LinearGroup,
                 *,
                 as_vector: "Optional[Vector]" = None,
                 index: Optional[int] = None,
                 mult: int = 1):
        self.G = G
        if as_vector is not None:
            assert len(as_vector) == G.rank
            self.as_vector = as_vector
        self.index = index
        self.mult = mult
    
    @cached_property
    def as_vector(self) -> "Vector":
        """ Vector representation of the weight """
        # Must be implemented in derived class if not given at construction.
        raise NotImplementedError

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Weight):
            return NotImplemented
        return self.G == other.G and self.as_vector == other.as_vector

    def leq(self,
            other: "Weight",
            symmetries: Optional[Iterable[int]] = None) -> bool:
        """
        Dominance order (partial sums) on each block given by sym. 

        Used in type boson and fermion
        The convention is 
        
            chi2.leq(chi,sym) and tau dominant => <chi2,tau> <= <chi2,tau>

        sym precise the sense of tau strictly dominant (decreasing on each
        block of sizes given by sym).
        """
        if symmetries is None:
            symmetries = self.G

        v = other.as_vector - self.as_vector
        shift = 0
        for s in symmetries:
            partial_sum = 0
            for i in range(s):
                partial_sum += v[shift + i]
                if partial_sum < 0:
                    return False
            if partial_sum != 0:
                return False
            shift += s
        return True

        """
        # FIXME: Alternative versions
        from .blocks import Blocks
        from itertools import accumulate
        diff = other.as_vector - self.as_vector
        for block in Blocks(diff, symmetries):
            for partial_sum in accumulate(block, initial=0):
                if partial_sum < 0:
                    return False
            if partial_sum != 0:
                return False
        return True
        """
            
    def __le__(self, other: "Weight") -> bool:
        """ Implementation of self <= other (partial ordering) """
        return self.leq(other)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(as_vector={self.as_vector}" + (f", idx: {self.index}" if self.index is not None else "") + ")"

    def orbit_symmetries(self, S: Iterable[int]) -> Iterable["Weight"]:
        return iter(())


class WeightAsList(Weight):
    """
    Weight class for the Kronecker representation

    Weight is stored as a list of one epsilon_i by copy of GL.
    """
    as_list: tuple[int, ...]

    def __init__(self, G: LinearGroup, as_list: Iterable[int], **kwargs: Any):
        super().__init__(G, **kwargs)
        self.as_list = tuple(as_list)
        assert (
            len(self.as_list) == len(self.G)
            and all(0 <= w < g for w, g in zip(self.as_list, self.G))
        ), "Invalid weight"

    @cached_property
    def as_vector(self) -> "Vector":
        from .rings import ZZ, vector
        v = vector(ZZ, self.G.rank)
        for shift, x in zip(itertools.accumulate(self.G, initial=0), self.as_list):
            v[shift + x] = 1
        return v

    def __le__(self, other: Weight) -> bool:
        """
        Implementation of self <= other (partial ordering)
        
        >>> from .representation import KroneckerRepresentation
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
        if not isinstance(other, WeightAsList):
            return NotImplemented
        return all(ws >= wo for ws, wo in zip(self, other))
    
    def leq(self,
            other: Weight,
            symmetries: Optional[Iterable[int]] = None) -> bool:
        if symmetries is None:
            return self <= other
        else:
            return super().leq(other, symmetries)
        
    def __eq__ (self, other: object) -> bool:
        if not isinstance(other, WeightAsList):
            return NotImplemented
        return self.G == other.G and self.as_list == other.as_list

    def __repr__(self) -> str:
        return f"WeightAsList({self.as_list}" + (f", idx: {self.index}" if self.index is not None else "") + ")"

    def __len__(self) -> int:
        return len(self.as_list)
    
    def __iter__(self) -> Iterator[int]:
        return iter(self.as_list)

    def orbit_symmetries(self, S: Iterable[int]) -> Iterable["Weight"]:
        from .utils import orbit_symmetries
        orbit = orbit_symmetries(self.as_list, S)
        for sym_w in orbit:
            yield WeightAsList(self.G, as_list=sym_w)


class WeightAsListOfList(Weight):
    """
    Weight class for the particle representation

    Weight is stored as a list of epsilon_i for each copy of GL.
    """
    as_list_of_list: tuple[tuple[int, ...], ...]

    def __init__(self, G: LinearGroup, as_list_of_list: Iterable[Iterable[int]], **kwargs: Any):
        super().__init__(G, **kwargs)
        self.as_list_of_list = tuple(tuple(w) for w in as_list_of_list)
        assert(
            len(self.as_list_of_list) == len(self.G)
            and all(all(0 <= wi < g for wi in w) for w, g in zip(self.as_list_of_list, self.G))
        ), "Invalid weight"
    
    @cached_property
    def as_vector(self) -> "Vector":
        from .rings import ZZ, vector
        v = vector(ZZ, self.G.rank)
        shift=0
        for shift, x in zip(itertools.accumulate(self.G, initial=0), self.as_list_of_list):
            for xi in x:
                v[shift + xi] += 1
        return v
    
    def __le__(self, other: Weight) -> bool:
        if not isinstance(other, WeightAsListOfList):
            return NotImplemented

        return all(
            all(ws >= wo for ws, wo in zip(ls, lo))
            for ls, lo in zip(self.as_list_of_list, other.as_list_of_list)
        )

    def leq(self,
            other: "Weight",
            symmetries: Optional[Iterable[int]] = None) -> bool:
        if symmetries is None:
            return self <= other
        else:
            return super().leq(other, symmetries)

    def __eq__ (self, other: object) -> bool:
        if not isinstance(other, WeightAsListOfList):
            return NotImplemented
        return self.G == other.G and self.as_list_of_list == other.as_list_of_list

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.as_list_of_list}" + (f", idx: {self.index}" if self.index is not None else "") + ")"