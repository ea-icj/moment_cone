from functools import cached_property

from .utils import prod, symmetries
from .typing import *


# Ugly things due to the circular import between
# Dimension and PolynomialRingForWeight
# FIXME: should be fixed if we separate Dimension from
# the associated polynomial rings, like in a kind of space class.
if TYPE_CHECKING:
    from .rings import PolynomialRingForWeights, Ring

__all__ = (
    "Dimension",
)

class Dimension(tuple[int, ...]):
    """
    Dimensions of the ??? space
    
    It now contains also some useful rings used in some parts of the module.

    TODO: it may be clearer to move all this in the kind of space class

    Examples:
    >>> d = Dimension((4, 4, 3, 2))
    >>> d
    (4, 4, 3, 2)
    >>> d.sum
    13
    >>> d.dimV
    96
    >>> d.dimU
    20
    >>> d.symmetries
    (2, 1, 1)

    It should also be noted that this class ensure uniqueness of an instance
    for a given sequence of dimensions:
    >>> d2 = Dimension((4, 4, 3, 2))
    >>> d is d2
    True
    """
    all_instances: dict["Dimension", "Dimension"] = {}

    def __new__(cls, dimensions):
        """ Construction with reusing of already computed Dimension instance """
        d = super().__new__(cls, dimensions)
        return cls.all_instances.setdefault(d, d)
    
    @cached_property
    def symmetries(self) -> tuple[int, ...]:
        """ Returns length of the symmetries in the dimensions """
        return tuple(symmetries(self))

    @cached_property
    def sum(self) -> int:
        return sum(self)

    @cached_property
    def dimV(self) -> int:
        """ Dimension of the vector space V """
        return prod(self)

    @cached_property
    def dimU(self) -> int:
        """ Dimension of the unipotent subgroup U """
        g = sum(i**2 for i in self)
        return (g - len(self)) // 2
    
    @cached_property
    def Q(self) -> "Ring":
        from .rings import QQ
        return QQ
    
    @cached_property
    def QI(self) -> "Ring":
        from .rings import QQ, I
        return QQ[I]

    @cached_property
    def QZ(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        return PolynomialRingForWeights(QQ, "z")
    
    @cached_property
    def QV(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        from .weight import Weight
        return PolynomialRingForWeights(QQ, weights=Weight.all(self))
    
    @cached_property
    def QV2(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ
        from .weight import Weight
        return PolynomialRingForWeights(
            QQ,
            weights=Weight.all(self),
            seed=('va', 'vb'),
        )

    @cached_property
    def QIV(self) -> "PolynomialRingForWeights":
        from .rings import PolynomialRingForWeights, QQ, I
        from .weight import Weight
        return PolynomialRingForWeights(
            QQ[I], 
            weights=Weight.all(self),
            seed=('vr', 'vi')
        )

    def u_max(self, e: "Dimension") -> int:
        """
        Maximal value of u obtained by extending a e-1-PS to a d-1-PS

        For a dimension vector d=(d_i), and a list of number of Levi blocks in each d_i,
        computes the maximal dimension of a nilradical.

        Examples:
        >>> d1 = Dimension((2, 3, 4))
        >>> d1.u_max(d1)
        10
        >>> d2 = Dimension((3, 5, 2))
        >>> d2.u_max(d2)
        14
        >>> d1.u_max(d2)
        8
        >>> d2.u_max(d1)
        11
        """
        from math import floor
        return sum(floor(di * di / 2 * (1 - 1 / ei)) for di, ei in zip(self, e))
        
        
        

