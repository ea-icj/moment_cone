"""
Tools to compute the dimension of the cone ???
"""

from sage.all import vector, I # type: ignore
from sage.structure.element import Vector # type: ignore

from .typing import *
from .weight import Weight
from .dimension import Dimension
from .polynomial_ring import PolynomialRingForWeights, Ring

def point_vect_QI(pds: Iterable[Weight], d: Dimension, ring: Ring, bounds: tuple[int, int] = (-100, 100)) -> Vector:
    """
    Generates a random vector of Vect(pds) over given ring.

    Coefficients are integers withing the given closed interval.
    """
    from random import randint
    v = vector(ring, d.prod)

    for chi in pds:
        coeffs = tuple(randint(*bounds) for _ in range(ring.degree()))
        if len(coeffs) == 1:
            value = ring(coeffs[0])
        else:
            value = ring(coeffs)
        v[chi.index_in(d)] = value

    return v

def point_vect_QV(pds: Iterable[Weight], d: Dimension, ring: PolynomialRingForWeights) -> Vector:
    """
    Generates a base vector of Vect(pds) over given polynomial (real) ring using the variables associated to the weights.
    """
    from .polynomial_ring import variable
    v = vector(ring.sage_ring, d.prod)

    for chi in pds:
        v[chi.index_in(d)] = ring.variable(chi)

    return v


def point_vect_QIV(pds: Iterable[Weight], d: Dimension, ring: PolynomialRingForWeights) -> Vector:
    """
    Generates a base vector of Vect(pds) over given polynomial (complex) ring using the variables associated to the weights.
    """
    from .polynomial_ring import variable
    v = vector(ring.sage_ring, d.prod)

    for chi in pds:
        vr, vi = ring.variable(chi)
        v[chi.index_in(d)] = vr + I * vi

    return v

def point_vect_QZ(pds: Iterable[Weight], d: Dimension, ring: Ring, bounds: tuple[int, int] = (-100, 100)) -> Vector:
    """
    Generates a random 1st-order polynomial of Vect(pds) over given ring.

    Coefficients are integers withing the given closed interval.
    """
    from random import randint
    from .polynomial_ring import variable
    v = vector(ring, d.prod)
    z = ring.variable("z")

    for chi in pds:
        coeffs = [randint(*bounds) for _ in range(2)]
        v[chi.index_in(d)] = coeffs[0] + z * coeffs[1]

    return v

def point_vect(pds: Iterable[Weight], d: Dimension, ring: PolynomialRingForWeights, bounds: tuple[int, int] = (-100, 100)) -> Vector:
    """ Returns an element of Vect(pds) following a method that depends on the given ring """
    match ring:
        case d.Q | d.QI:
            return point_vect_QI(pds, d, ring, bounds)
        case d.QZ:
            return point_vect_QZ(pds, d, ring, bounds)
        case d.QV:
            return point_vect_QV(pds, d, ring)
        case d.QIV:
            return point_vect_QIV(pds, d, ring)
        case _:
            raise ValueError("Unknown ring")
        
def rank_RC(M, d: Dimension) -> int:
    """
    Rank of the R-linear from R^(number of columns) to C^(number of rows)
    """
    return NotImplemented