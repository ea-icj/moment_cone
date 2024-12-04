from functools import cached_property

from .utils import prod
from .typing import *

# Ugly things due to the circular import between
# Dimension and PolynomialRingForWeight
# FIXME: should be fixed if we separate Dimension from
# the associated polynomial rings, like in a kind of space class.
if TYPE_CHECKING:
    from .polynomial_ring import PolynomialRingForWeights, Ring

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
    >>> d.prod
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
    def prod(self) -> int:
        return prod(self)
    
    @cached_property
    def Q(self) -> "Ring":
        from sage.all import QQ # type: ignore
        return QQ
    
    @cached_property
    def QI(self) -> "Ring":
        from sage.all import QQ, I # type: ignore
        return QQ[I]

    @cached_property
    def QZ(self) -> "PolynomialRingForWeights":
        from .polynomial_ring import PolynomialRingForWeights
        from sage.all import QQ # type: ignore
        return PolynomialRingForWeights(QQ, "z")
    
    @cached_property
    def QV(self) -> "PolynomialRingForWeights":
        from .polynomial_ring import PolynomialRingForWeights
        from .weight import Weight
        from sage.all import QQ # type: ignore
        return PolynomialRingForWeights(QQ, weights=Weight.all(self))
    
    @cached_property
    def QV2(self) -> "PolynomialRingForWeights":
        from .polynomial_ring import PolynomialRingForWeights
        from .weight import Weight
        from sage.all import QQ # type: ignore
        return PolynomialRingForWeights(
            QQ,
            weights=Weight.all(self),
            seed=('va', 'vb'),
        )

    @cached_property
    def QIV(self) -> "PolynomialRingForWeights":
        from .polynomial_ring import PolynomialRingForWeights
        from .weight import Weight
        from sage.all import QQ, I # type: ignore
        return PolynomialRingForWeights(
            QQ[I], 
            weights=Weight.all(self),
            seed=('vr', 'vi')
        )
