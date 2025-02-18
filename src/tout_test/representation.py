__all__ = (
    "Representation",
    "KroneckerRepresentation",
    "ParticleRepresentation",
)

from abc import ABC, abstractmethod
from functools import cached_property
import itertools

from .typing import *
from .linear_group import LinearGroup
from .weight import Weight as WeightBase, WeightAsList, WeightAsListOfList


class Representation(ABC):
    """ Base class of a representation """
    Weight: ClassVar[type[WeightBase]] = WeightBase # Weight class
    G: LinearGroup

    def __init__(self, G: LinearGroup):
        self.G = G

    def weight(self, *args, **kwargs) -> WeightBase:
        """ Creates a weight for the given representation """
        return self.Weight(self.G, *args, **kwargs)
    
    @abstractmethod
    def all_weights(self) -> Iterable[WeightBase]:
        ...

    @abstractmethod
    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        ...

    @cached_property
    @abstractmethod
    def dim_cone(self) -> int:
        """ Expected dimension of the cone. To be checked with Stabilizer of K """
        ...

    @cached_property
    @abstractmethod
    def dim(self) -> int:
        """ Dimension of V """
        ...


class KroneckerRepresentation(Representation):
    Weight = WeightAsList
    
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
            yield WeightAsList(self.G, w, index=idx)

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        """
        Returns index of the weight in the lexicographical order for the current representation (see `all_weights` method).
        
        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """
        if not use_internal_index or chi.index is None:
            if not isinstance(chi, WeightAsList):
                raise ValueError("Invalid weight representation")        
            from operator import mul
            stride = itertools.accumulate(reversed(self.G[1:]), mul, initial=1)
            chi.index = sum(v * s for v, s in zip(reversed(chi.as_list), stride))
        return chi.index


class ParticleRepresentation(Representation):
    """ Representation specific to physical particles """
    particle_cnt: int

    def __init__(self, G: LinearGroup, particle_cnt: int):
        super().__init__(G)
        self.particle_cnt = particle_cnt
        if len(G) != 1:
            raise NotImplementedError("Product of GL not supported for particle representation")

    @cached_property
    def dim_cone(self) -> int:
        return self.G.rank
    

class FermionRepresentration(ParticleRepresentation):
    """ Representation specific to Fermion particles """
    Weight = WeightAsListOfList

    @cached_property
    def dim(self) -> int:
        from math import comb
        return comb(self.G.rank, self.particle_cnt)

    def all_weights(self) -> Iterable[WeightBase]:
        weights = itertools.combinations(range(self.G[0]), self.particle_cnt)
        for i, w in enumerate(weights):
            yield WeightAsListOfList(self.G, as_list_of_list=(w,), index=i)

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        if not isinstance(chi, WeightAsListOfList):
            raise ValueError("Invalid weight representation")
        
        if not use_internal_index or chi.index is None:
            raise NotImplementedError() # TODO
        
        return chi.index
    

class BosonRepresentration(ParticleRepresentation):
    """ Representation specific to Boson particles """
    Weight = WeightAsListOfList

    @cached_property
    def dim(self) -> int:
        from math import comb
        return comb(
            self.G.rank + self.particle_cnt - 1,
            self.G.rank - 1
        )