import itertools

from .typing import *
from .kronecker import KroneckerCoefficient
from .partition import Partition
from .permutation import Permutation

def all_partitions_of_max_length(n: int, l: Sequence[int], kro: KroneckerCoefficient) -> Iterable[tuple[tuple[Partition, ...], int]]:
    """
    All len(l)-uplets of partitions of n of non-zero Kronecker coefficient
    and so that len(p_i) <= l_i
    """
    # Sort by increasing maximal length (faster) and keep order
    permutation, sorted_l = zip(*sorted(enumerate(l), key=lambda l: l[1]))
    p_inverse = Permutation(permutation).inverse

    # All nuplet of partitions without the last length constraint
    head_product = itertools.product(
        *(Partition.all_for_integer(n, li) for li in sorted_l[:-1])
    )

    # Computing the product and yielding only the partitions of the decomposition
    # whose length respect the last constraint.
    for head in head_product:
        product = kro.product(head)
        for p, c in product.items():
            if len(p) <= sorted_l[-1]:
                yield p_inverse(head + (p,)), c
