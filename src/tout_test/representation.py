from abc import ABC, abstractmethod

from .typing import *
from .linear_group import LinearGroup
from .weight import Weight as WeightBase

class Representation(ABC):
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

    @property
    @abstractmethod
    def dim_cone(self) -> int:
        """ Expected dimension of the cone. To be checked with Stabilizer of K """
        ...

    @property
    @abstractmethod
    def dim(self) -> int:
        """ Dimension of V """
        ...





