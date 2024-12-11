from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .cone_dimension import point_vect, vector, matrix
from .spaces import action_op_el

# FIXME: we get d from tau but in the current code, it will leads to recreate the rings for each tau.
def is_not_contracted(inverse_v: Sequence[Root], tau: Tau, method: Method) -> bool:
    """
    ???

    The vector v is chosen depending on the method: probabilistic or symbolic.

    FIXME: example

    Example:
    >>> from cone import *
    >>> d = Dimension((2, 3, 4))
    >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
    >>> tau
    1 | 6 2 | 1 4 1 | 4 5 3 1
    >>> roots = list(Root.all_of_U(d))[:10]
    >>> is_not_contracted(roots, tau, "probabilistic")
    False
    >>> is_not_contracted(roots, tau, "symbolic")
    False
    """
    d = tau.d

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = d.QI
    elif method == "symbolic":
        ring = d.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    
    # FIXME: do we stay we list conversion at each call?
    # Maybe grading root and weight should be implemented using a more convenient class?
    from itertools import chain
    non_positive_weights = list(chain.from_iterable(tau.non_positive_weights.values()))
    positive_weights = list(chain.from_iterable(tau.positive_weights.values()))

    # FIXME: why not positive_weights instead of non_negative_weights?
    v = point_vect(non_positive_weights, d, ring, bounds=(-1000, 1000))

    T = matrix(ring, len(positive_weights), len(inverse_v))
    for j, root in enumerate(inverse_v):
        uv = action_op_el(root, v, d)
        for i, chi in enumerate(positive_weights):
            T[i, j] = uv[chi.index_in(d)]

    rank_T = T.change_ring(ring.fraction_field()).rank()
    return rank_T == len(inverse_v)


