from .typing import *
from .dimension import Dimension
from dataclasses import dataclass
import itertools

__all__ = (
    "Root",
)

@dataclass(frozen=True, slots=True)
class Root:
    """ Root element for tau """
    k: int
    i: int
    j: int

    @property
    def is_in_U(self) -> bool:
        """ Check if this root is in U """
        return self.i < self.j
    
    @staticmethod
    def all_of_U(d: Dimension) -> Iterable["Root"]:
        """ Returns all possible root from U for given dimensions """
        for k, dk in enumerate(d):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)

    @staticmethod
    def all(d: Dimension) -> Iterable["Root"]:
        """ Returns all possible roots from G for given dimensions """
        for k, dk in enumerate(d):
            for i, j in itertools.product(range(dk), repeat=2) :
                if i != j:
                    yield Root(k, i, j)

    @staticmethod
    def all_of_T(d: Dimension) -> Iterable["Root"]:
        """ Returns all possible roots from T for given dimensions """
        for k, dk in enumerate(d):
            for i in range(dk):
                yield Root(k, i, i)            

