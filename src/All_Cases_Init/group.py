from functools import cached_property
#from typing import Generic,TypeVar
#T=TypeVar('T')
#from collections.abc import Iterable, Sequence, Callable

from .utils import prod, symmetries
from .typing import *


__all__ = (
    "LinGroup",
)

class LinGroup(tuple[int, ...]):
    """
    Product of Linear Groups GL(d_i)

    Examples:
    >>> G = LinGroup([4, 4, 3, 2])
    >>> G
    GL(4)xGL(4)xGL(3)xGL(2)
    >>> G.rank
    13
    >>> G.dim
    45
    >>> G.dimU
    16
    >>> G.outer
    (2,1,1)
    It should also be noted that this class ensure uniqueness of an instance
    for a given sequence of dimensions:
    >>> G2=LinGroup([4,4,3,2])
    >>> G==G2
    True
    """
    all_instances: dict["LinGroup", "LinGroup"] = {}

    def __new__(cls, dimensions):
        """ Construction with reusing of already computed Dimension instance """
        d = super().__new__(cls, dimensions)
        return cls.all_instances.setdefault(d, d)

    def __repr__(self) -> str:
        return 'x'.join(f'GL({i})' for i in self)

    @cached_property
    def outer(self) -> tuple[int, ...]:
        """ Returns length of the symmetries in the dimensions """
        return tuple(symmetries(self))

    @cached_property
    def rank(self) -> int:
        """
        The rank of G
        """
        return sum(self)

    @cached_property
    def dim(self) -> int:
        """ 
        Dimension of the group G 
        """
        return sum(i**2 for i in self)

    @cached_property
    def dimU(self) -> int:
        """ 
        Dimension of the unipotent subgroup U 
        """
        g = sum(i**2 for i in self)
        return (self.dim - self.rank) // 2

    def u_max(self, Gred: "LinGroup") -> int:
        """
        Maximal value of u obtained by extending a e-1-PS to a d-1-PS

        For a Linear Group d=(d_i), and a list of number of Levi blocks in each e_i,
        computes the maximal dimension of a nilradical.
        
        Maximal is relative to the various embeddings of L in self
        
        Warning : need d_i and e_i given ordered
        Examples:
        >>> G1 = LinGroup([5,3,2])
        >>> G1.u_max(G1)
        10
        >>> G2 = LinGroup([2, 2, 1])
        >>> G2.u_max(G2)
        2
        >>> G1.u_max(G2)
        8
        >>> G1.u_max(G1)
        14
        """
        from math import floor
        return sum(floor(d * d / 2 * (1 - 1 / e)) for d, e in zip(self, Gred))
    
    def QU(self, base_ring=None) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        from .root import Root
        if base_ring==None:
            base_ring=QQ
        return PolynomialRingForWeights(base_ring,
            weights=Root.all_of_U(self),
            seed=('u'),
        )
    
    @cached_property
    def QU_Q(self) -> "PolynomialRingForWeights":
        return self.QU()
        
    

        

