"""
Tools to compute the dimension of the cone ???
"""

from .typing import *
from .weight import Weight
from .dimension import Dimension
from .tau import Tau
from .root import Root
from .rings import PolynomialRingForWeights, Ring, vector, Vector, matrix, Matrix, I
from .rings import real_part, imag_part
from .spaces import action_op_el

def point_vect_QI(pds: Iterable[Weight], d: Dimension, ring: Ring, bounds: tuple[int, int] = (-100, 100)) -> Vector:
    """
    Generates a random vector of Vect(pds) over given ring.

    Coefficients are integers withing the given closed interval.
    """
    from random import randint
    v = vector(ring, d.dimV)

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
    v = vector(ring, d.dimV)

    for chi in pds:
        v[chi.index_in(d)] = ring.variable(chi)

    return v

def point_vect_QV2(pds: Iterable[Weight], d: Dimension, ring: PolynomialRingForWeights) -> Vector:
    """
    Generates a base vector of Vect(pds) over given polynomial (real) ring using the variables associated to the weights.

    For each weights, at its index, returns va_chi * z + vb_chi
    """
    v = vector(ring, d.dimV)
    z = ring.variable("z")
    for chi in pds:
        va, vb = ring.variable(chi)
        v[chi.index_in(d)] = va * z + vb
    return v

def point_vect_QIV(pds: Iterable[Weight], d: Dimension, ring: PolynomialRingForWeights) -> Vector:
    """
    Generates a base vector of Vect(pds) over given polynomial (complex) ring using the variables associated to the weights.
    """
    v = vector(ring, d.dimV)

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
    v = vector(ring, d.dimV)
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
        case d.QV2:
            return point_vect_QV2(pds, d, ring)
        case d.QIV:
            return point_vect_QIV(pds, d, ring)
        case _:
            raise ValueError("Unknown ring")
        
def rank_RC(M: Matrix) -> int:
    """
    Rank of the R-linear from R^(number of columns) to C^(number of rows)

    >>> d = Dimension((4, 3, 2))
    >>> weights = list(Weight.all(d))[5:12]
    >>> M = matrix(d.QIV, 10, 10)
    >>> vr, vi = d.QIV.variable(weights[2])
    >>> M[0, 0] = vr + I * vi
    >>> M[0, 5] = 3 * vr
    >>> rank_RC(M)
    2
    >>> M.change_ring(M.base_ring().fraction_field()).rank()
    1
    """
    N = matrix(M.base_ring(), 2 * M.nrows(), M.ncols())
    N[:M.nrows(), :] = real_part(M)
    N[M.nrows():, :] = imag_part(M)

    return N.change_ring(M.base_ring().fraction_field()).rank()


# FIXME: we get d from tau but in the current code, it will leads to recreate the rings for each tau.
# TODO: verify example (and maybe add one better)
def dim_of_stabilizer_in_K_tau(tau: Tau, method: Method) -> int:
    """
    Compute the dimension of a stabilizer in $K^\tau$ for a vector v.
    
    The vector v is chosen depending on the method: probabilistic or symbolic.

    Here, K is the product of the GL.

    Example:
    >>> d = Dimension((2, 3, 4))
    >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
    >>> dim_of_stabilizer_in_K_tau(tau, "probabilistic")
    11
    >>> dim_of_stabilizer_in_K_tau(tau, "symbolic")
    11
    """
    d = tau.d # FIXME: Here, d is recreated from scratch, without rings. Should we ensure the uniqueness of the instance of d?

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = d.QI
    elif method == "symbolic":
        ring = d.QIV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    
    roots_K_tau = tau.orthogonal_roots
    weights_K_tau = tau.orthogonal_weights

    v = point_vect(weights_K_tau, d, ring, bounds=(-1000, 1000))

    dim_K_tau = d.sum + 2 * len(roots_K_tau)
    M = matrix(ring, len(weights_K_tau), dim_K_tau)
    j = 0

    # First block
    for k, dk in enumerate(d):
        for l in range(dk):
            root = Root(k, l, l) # FIXME: create a Root.all_diagonal?
            tv = action_op_el(root, v, d)
            for i, chi in enumerate(weights_K_tau):
                M[i, j] = tv[chi.index_in(d)]
            j += 1

    # Second block
    for root in roots_K_tau:
        tv_pos = action_op_el(root, v, d)
        tv_neg = action_op_el(Root(root.k, root.j, root.i), v, d) # FIXME: root.transpose?
        tv_k1 = tv_pos - tv_neg
        tv_k2 = I * (tv_pos + tv_neg)
        for i, chi in enumerate(weights_K_tau):
            M[i, j] = tv_k1[chi.index_in(d)]
            M[i, j + 1] = tv_k2[chi.index_in(d)]
        j += 2

    return dim_K_tau - rank_RC(M)

    

