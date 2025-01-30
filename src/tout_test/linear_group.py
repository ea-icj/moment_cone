__all__ = (
    "LinearGroup",
)

from functools import cached_property

from .utils import symmetries
from .typing import *

class LinearGroup(tuple[int, ...]):
    """
    Product of Linear Groups GL(d_i)

    Examples:
    >>> G = LinearGroup((4, 4, 3, 2))
    >>> G
    GL(4)xGL(4)xGL(3)xGL(2)
    >>> G.rank
    13
    >>> G.dim
    45
    >>> G.dimU # FIXME: expected 20 ?! but got 16...
    16
    >>> G.outer # FIXME: expected (1, 3) but got (2, 1, 1)
    (2, 1, 1)

    It should also be noted that this class ensure uniqueness of an instance
    for a given sequence of dimensions:
    >>> G2 = LinearGroup((4, 4, 3, 2))
    >>> G == G2
    True
    """
    all_instances: ClassVar[dict["LinearGroup", "LinearGroup"]] = {}

    def __new__(cls, dimensions: Iterable[int] ):
        """ Construction with reusing of already computed LinearGroupe instance """
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
        return sum(self)

    @cached_property
    def dim(self) -> int:
        """ Rank of the group G """
        return sum(i**2 for i in self)

    @cached_property
    def dimU(self) -> int:
        """ Dimension of the unipotent subgroup U """
        g = sum(i**2 for i in self)
        return (self.dim - self.rank) // 2

        

