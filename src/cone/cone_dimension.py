from sage.all import vector, Ring, RingElement # type: ignore
from sage.structure.element import Vector # type: ignore

from random import randint

from .typing import *
from .weight import Weight
from .dimension import Dimension

def PointV(pds: Iterable[Weight], d: Dimension, ring: Ring, upper: int = 1000) -> Vector:
    """ Generate a random vector of Vect(pds) over given ring """
    v = vector(ring, d.prod)

    for chi in pds:
        coeffs = tuple(randint(1, upper) for _ in range(ring.degree()))
        if len(coeffs) == 1:
            value = ring(coeffs[0])
        else:
            value = ring(coeffs)
        v[chi.index_in(d)] = value

    return v