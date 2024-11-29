"""
Tools to compute the dimension of the cone ???
"""

from sage.all import vector, Ring, RingElement, I # type: ignore
from sage.structure.element import Vector # type: ignore

from .typing import *
from .weight import Weight
from .dimension import Dimension

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

def point_vect_QV(pds: Iterable[Weight], d: Dimension, ring: Ring) -> Vector:
    """
    Generates a base vector of Vect(pds) over given polynomial (real) ring using the variables associated to the weights.
    """
    from .polynomial_ring import variable
    v = vector(ring, d.prod)

    for chi in pds:
        v[chi.index_in(d)] = variable(ring, chi)

    return v


def point_vect_QIV(pds: Iterable[Weight], d: Dimension, ring: Ring) -> Vector:
    """
    Generates a base vector of Vect(pds) over given polynomial (complex) ring using the variables associated to the weights.
    """
    from .polynomial_ring import variable
    v = vector(ring, d.prod)

    for chi in pds:
        v[chi.index_in(d)] = variable(ring, chi, seed="vr") + I * variable(ring, chi, seed="vi")

    return v

def point_vect_QZ(pds: Iterable[Weight], d: Dimension, ring: Ring, bounds: tuple[int, int] = (-100, 100)) -> Vector:
    """
    Generates a random 1st-order polynomial of Vect(pds) over given ring.

    Coefficients are integers withing the given closed interval.
    """
    from random import randint
    from .polynomial_ring import variable
    v = vector(ring, d.prod)
    z = variable(ring, "z")

    for chi in pds:
        coeffs = [randint(*bounds) for _ in range(2)]
        v[chi.index_in(d)] = coeffs[0] + z * coeffs[1]

    return v

def point_vect(pds: Iterable[Weight], d: Dimension, ring: Ring, bounds: tuple[int, int] = (-100, 100)) -> Vector:
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