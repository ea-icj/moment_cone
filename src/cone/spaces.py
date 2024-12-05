"""
Space related functions

Currently only for Kronecker. Later for Fermions and Bosons
"""

from .typing import *
from .root import Root
from .dimension import Dimension
from .rings import vector, Vector
from .weight import Weight

def action_op_el(alpha: Root, v: Vector, d: Dimension) -> Vector:
    """
    Elementary action of a root on the vector v
    
    Examples:
    >>> from cone import *
    >>> d = Dimension((2, 3))
    >>> v = vector(QQ, d.dimV)
    >>> v[:] = range(1, d.dimV + 1)
    >>> v
    (1, 2, 3, 4, 5, 6)
    >>> alpha = Root(k=1, i=0, j=2)
    >>> action_op_el(alpha, v, d)
    (3, 0, 0, 6, 0, 0)

    >>> d = Dimension((2, 3, 4))
    >>> v = vector(QQ, d.dimV)
    >>> v[:] = range(1, d.dimV + 1)
    >>> alpha = Root(k=1, i=0, j=2)
    >>> action_op_el(alpha, v, d)
    (9, 10, 11, 12, 0, 0, 0, 0, 0, 0, 0, 0, 21, 22, 23, 24, 0, 0, 0, 0, 0, 0, 0, 0)
    """
    assert len(v) == d.dimV
    vp = vector(v.base_ring(), d.dimV)

    d_not_k = list(d)
    del d_not_k[alpha.k]

    for chi in Weight.all(Dimension(d_not_k)):
        chi_i = list(chi)
        chi_i.insert(alpha.k, alpha.i)
        chi_j = list(chi)
        chi_j.insert(alpha.k, alpha.j)
        vp[Weight(chi_i).index_in(d)] = v[Weight(chi_j).index_in(d)]

    return vp




