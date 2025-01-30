__all__ = (
    'ParticleWeight',
    'ParticleRepresentation',
)

from functools import cached_property
import itertools

from .typing import *
from .representation import Representation
from .weight import Weight as WeightBase
from .linear_group import LinearGroup
from .rings import Vector, vector

class ParticleWeight(WeightBase):
    """
    Weight class for the particle representation
    """
    __weights: tuple[tuple[int, ...], ...]

    def __init__(self, G: LinearGroup, weights: Iterable[Iterable[int]], **kwargs):
        super().__init__(G, **kwargs)
        self.__weights = tuple(tuple(w) for w in weights)
        assert(
            len(self.__weights) == len(self.G)
            and all(all(0 <= wi < g for wi in w) for w, g in zip(self.__weights, self.G))
        ), "Invalid weight"
    
    @cached_property
    def as_vector(self) -> Vector:
        from .rings import ZZ
        v = vector(ZZ, self.G.rank)
        shift=0
        for shift, x in zip(itertools.accumulate(self.G, initial=0), self.__weights):
            for xi in x:
                v[shift + xi] += 1
        return v
    
    def __le__(self, other: WeightBase) -> bool:
        if not isinstance(other, ParticleWeight):
            return NotImplemented

        return all(
            all(ws >= wo for ws, wo in zip(ls, lo))
            for ls, lo in zip(self.__weights, other.__weights)
        )

    def leq(self,
            other: "WeightBase",
            symmetries: Optional[Iterable[int]] = None) -> bool:
        if symmetries is None:
            return self <= other
        else:
            return super().leq(other, symmetries)

    def __eq__ (self, other: object) -> bool:
        if not isinstance(other, ParticleWeight):
            return NotImplemented
        return self.G == other.G and self.__weights == other.__weights


    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.__weights}" + (f", idx: {self.index}" if self.index is not None else "") + ")"


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
