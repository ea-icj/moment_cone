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
    
    @property
    def opposite(self) -> "Root":
        """ Return the opposite of a root (k,i,j -> kj,i) """
        return Root(self.k, self.j, self.i)
    
    @staticmethod
    def all_of_U(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible root from U for given dimensions
        
        Example:
        >>> d = Dimension((2, 3))
        >>> for root in Root.all_of_U(d):
        ...     print(root)
        Root(k=0, i=0, j=1)
        Root(k=1, i=0, j=1)
        Root(k=1, i=0, j=2)
        Root(k=1, i=1, j=2)
        """
        for k, dk in enumerate(d):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)

    @staticmethod
    def all(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible roots from G (i != j) for given dimensions
        
        Example:
        >>> d = Dimension((2, 3))
        >>> for root in Root.all(d):
        ...     print(root)
        Root(k=0, i=0, j=1)
        Root(k=0, i=1, j=0)
        Root(k=1, i=0, j=1)
        Root(k=1, i=1, j=0)
        Root(k=1, i=0, j=2)
        Root(k=1, i=2, j=0)
        Root(k=1, i=1, j=2)
        Root(k=1, i=2, j=1)
        """
        for k, dk in enumerate(d):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)
                yield Root(k, j, i)


    @staticmethod
    def all_of_T(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible roots from T (i == j) for given dimensions
        
        
        Example:
        >>> d = Dimension((2, 3))
        >>> for root in Root.all_of_T(d):
        ...     print(root)
        Root(k=0, i=0, j=0)
        Root(k=0, i=1, j=1)
        Root(k=1, i=0, j=0)
        Root(k=1, i=1, j=1)
        Root(k=1, i=2, j=2)
        """
        for k, dk in enumerate(d):
            for i in range(dk):
                yield Root(k, i, i)            

