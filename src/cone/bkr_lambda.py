import itertools

from numpy.typing import NDArray

from .typing import *
from .kronecker import KroneckerCoefficient
from .partition import Partition
from .permutation import Permutation
from .utils import prod


def all_partitions_of_max_length(n: int, l: Sequence[int], kro: KroneckerCoefficient) -> Iterable[tuple[tuple[Partition, ...], int]]:
    """
    All len(l)-uplets of partitions of n of non-zero Kronecker coefficient
    and so that len(p_i) <= l_i

    Examples:
    >>> from cone.kronecker import KroneckerCoefficientMLCache
    >>> kc = KroneckerCoefficientMLCache()
    >>> for p in all_partitions_of_max_length(3, [2, 3, 2], kc):
    ...    print(p)
    ((Partition((3,)), Partition((3,)), Partition((3,))), 1)
    ((Partition((3,)), Partition((2, 1)), Partition((2, 1))), 1)
    ((Partition((2, 1)), Partition((2, 1)), Partition((3,))), 1)
    ((Partition((2, 1)), Partition((1, 1, 1)), Partition((2, 1))), 1)
    ((Partition((2, 1)), Partition((2, 1)), Partition((2, 1))), 1)
    ((Partition((2, 1)), Partition((3,)), Partition((2, 1))), 1)

    >>> for p in all_partitions_of_max_length(3, [2, 2], kc):
    ...    print(p)
    ((Partition((3,)), Partition((3,))), 1)
    ((Partition((2, 1)), Partition((2, 1))), 1)
    >>> list(all_partitions_of_max_length(3, [2, 0, 2], kc))
    []
    >>> list(all_partitions_of_max_length(3, [2], kc))
    [((Partition((3,)),), 1)]
    >>> list(all_partitions_of_max_length(3, [0], kc))
    []
    >>> list(all_partitions_of_max_length(3, [], kc))
    []
    """
    if len(l) == 1 and l[0] >= 1:
        yield (Partition(n),), 1
    elif len(l) == 2:
        for p in Partition.all_for_integer(n, max_length=min(l)):
            yield (p, p), 1
    elif len(l) >= 3:
        # Sort by increasing maximal length (faster) and keep order
        permutation, sorted_l = zip(*sorted(enumerate(l), key=lambda l: l[1]))
        p_inverse = Permutation(permutation).inverse

        # All nuplets of partitions without the last length constraint
        head_product = itertools.product(
            *(Partition.all_for_integer(n, li) for li in sorted_l[:-1])
        )

        # Computing the product and yielding only the partitions of the decomposition
        # whose length respects the last constraint.
        for head in head_product:
            product = kro.product(head)
            for p, c in product.items():
                if len(p) <= sorted_l[-1]:
                    yield p_inverse(head + (p,)), c


def all_lambda_matrix(delta: Sequence[int], max_length: NDArray, kro: KroneckerCoefficient) -> Iterable[tuple[NDArray, int]]:
    """ All Lambda matrices form given weight vector and maximal length constraints
    
    Yield a matrix and the product of the Kronecker coefficient of each row.
    """
    N, s = max_length.shape
    assert N == len(delta)

    import numpy as np

    row_product = itertools.product(*(
        all_partitions_of_max_length(n, l, kro)
        for n, l in zip(delta, max_length)
    ))

    for rows_and_coeff in row_product:
        lambda_matrix = np.empty((N, s), dtype=object)
        lambda_matrix[:, :], coeffs = zip(*rows_and_coeff)
        yield lambda_matrix, prod(coeffs)

