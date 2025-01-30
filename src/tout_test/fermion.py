from functools import cached_property
import itertools

from .typing import *
from .particle import ParticleRepresentation, ParticleWeight
from .weight import Weight as WeightBase
from .linear_group import LinearGroup

class FermionWeight(ParticleWeight):
    """
    Weight class for the Fermion representation
    """



class FermionRepresentration(ParticleRepresentation):
    Weight = FermionWeight

    @cached_property
    def dim(self) -> int:
        from math import comb
        return comb(self.G.rank, self.particle_cnt)

    def all_weights(self) -> Iterable[WeightBase]:
        weights = itertools.combinations(range(self.G[0]), self.particle_cnt)
        for i, w in enumerate(weights):
            yield FermionWeight(self.G, weights=(w,), index=i)

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        if not isinstance(chi, FermionWeight):
            raise ValueError("Invalid weight representation")
        
        if not use_internal_index or chi.index is None:
            raise NotImplementedError() # TODO
        
        return chi.index