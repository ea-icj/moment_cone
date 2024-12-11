from functools import cached_property

from .utils import prod
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
    >>> d.symmetries
    (2, 1, 1)
    """
    @cached_property
    def symmetries(self) -> tuple[int, ...]:
        """ Returns length of the symmetries in the dimensions """
        from .utils import group_by_block
        return tuple(length for _, length in group_by_block(self))

    @cached_property
    def sum(self) -> int:
        return sum(self)

    @cached_property
    def dimV(self) -> int:
        """ Dimension of the vectorial space V """
        return prod(self)
    
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
            "z",
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

    def uMAX(self,e: "Dimension")->int: # Maximal value of u obtained by extending a e-1-PS to a d-1-PS
        """
        For a dimension vector d=(d_i), and a list of number of Levi blocks in each d_i,
        computes the maximal dimension of a nilradical.
        """
        from math import floor
        return(sum([floor(self[i]*self[i]/2*(1-1/e[i])) for i in range(len(self))]))
        
        

