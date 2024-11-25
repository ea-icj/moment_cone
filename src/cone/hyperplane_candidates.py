from dataclasses import dataclass

from .typing import *
from .weight import Weight
from .dimension import Dimension
from .utils import short_prod

from sage.all import matrix, ZZ

__all__ = (
    "find_hyperplanes",
)


@dataclass(slots=True)
class WeightSieve:
    """ Sieve used to categorized weights during the search of hyperplane """
    indeterminate: list[Weight] # Weights whose status is not currently determinate
    excluded: list[Weight] # Weights excluded from being zero
    zero: list[Weight] # Weights considered as zero
    positive: list[Weight] # Positive weights
    negative: list[Weight] # Negative weights

    def copy(self) -> "WeightSieve":
        return WeightSieve(
            self.indeterminate.copy(),
            self.excluded.copy(),
            self.zero.copy(),
            self.positive.copy(),
            self.negative.copy(),
        )


def smart_remove(l: list[T], idx: int):
    """
    Remove an element of a list by swapping it with the last element

    It saves the O(N) cost of shift all the remaining elements.

    TODO: if necessary, it could be possible to use this trick
    to remove the necessary of copying each list by always moving removed
    elements to the end and by emulating a shorter list while being able
    to recover the previous size (and thus the removed elements).
    """
    if idx < len(l) - 1:
        l[idx] = l.pop()
    else:
        l.pop()

def sign_assignment(chi: Weight,
                    S_input: list[Weight],
                    S_negative: list[Weight],
                    S_positive: list[Weight]):
    """ Determining the sign of weight from S_input by comparing it to the curring weight chi """
    idx = 0
    while idx < len(S_input):
        chi2 = S_input[idx]
        if chi2 <= chi:
            S_negative.append(chi2)
            smart_remove(S_input, idx)
        elif chi <= chi2:
            S_positive.append(chi2)
            smart_remove(S_input, idx)
        else:
            idx += 1 # index incremented only when not element where removed

def hyperplane_matrix(S: Sequence[Weight], d: Dimension):
    """ Matrix with columns the weights indexed by the set S """
    M = matrix(ZZ, d.sum + 1, len(S))
    for j, chi in enumerate(S):
        M[0, j] = 1
        shift = 1
        for u in range(len(chi)):        
            M[shift + chi[u], j] = 1
            shift += d[u]
    return M

def check_hyperplane_dim(S: Sequence[Weight], d: Dimension):
    """ Check that the given sequence of weight has the dimension of an hyperplane """
    target_rank = d.sum - len(d)
    if len(S) < target_rank:
        return False
    
    M = hyperplane_matrix(S, d)
    return M.rank(algorithm="flint") == target_rank # Flint algorithm is faster than the default one

def has_enough_leq_weights(chi: Weight, u: int) -> bool:
    """ True for weights with more that u weights bigger for the order leq """
    leq_cnt = short_prod(c + 1 for c in chi) - 1
    # FIXME: verify this function, it differ from the original version
    return leq_cnt > u

def simplify_regularity(St: WeightSieve, d: Dimension, chi: Weight) -> None:
    """
    Changing a single epsilon^s_i into epsilon^s_j prevents regularity of tau, since this introduces equation epsilon^s_i-epsilon^s_j. When c gets in the zero list, the function removes these weights (with respect to c) from the indeterminate liste and places them in the zero list.
    """
    for k, dk in enumerate(d):
        mchi2 = list(chi)
        for y in range(dk):
            if y != chi[k]:
                mchi2[k] = y
                chi2 = Weight(mchi2)
                try:
                    idx = St.indeterminate.index(chi2)
                except ValueError:
                    pass
                else:
                    smart_remove(St.indeterminate, idx)
                    St.excluded.append(chi2)

def find_hyperplanes(d: Dimension, u: int) -> Iterable[list[Weight]]:
    """
    Returns hyperplane candidates

    What is u?
    """
    St = WeightSieve([], [], [], [], [])
    for chi in Weight.all(d):
        if has_enough_leq_weights(chi, u):
            St.positive.append(chi)
        else:
            St.indeterminate.append(chi)
    return find_hyperplanes_impl(St, d, u)

def find_hyperplanes_impl(St: WeightSieve, d: Dimension, u: int) -> Iterable[list[Weight]]:
    """ Recursive part to find the hyperplane candidates """
    if check_hyperplane_dim(St.zero, d):
        # Candidate hyperplane if the dimension is appropriate
        yield St.zero

    elif len(St.zero) + len(St.indeterminate) >= d.sum - len(d) and len(St.indeterminate) > 0:
        # Next element to consider
        chi = St.indeterminate.pop()

        # Two possible actions with this element:

        # 1. We exclude it from the possible zero elements
        St.excluded.append(chi)
        yield from find_hyperplanes_impl(St, d, u)
        St.excluded.pop()

        # 2. We define it as a zero element (on the hyperplane)
        St2 = St.copy()
        St2.zero.append(chi)

        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(chi, St2.indeterminate, St2.negative, St2.positive)
        sign_assignment(chi, St2.excluded, St2.negative, St2.positive)

        # 2.2 Continuing if there are not too much negative elements
        if len(St2.negative) <= u:
            # simplify_regularity(St2, d, chi) # Using regularity condition of tau
            # TODO: check if it is necessary
            yield from find_hyperplanes_impl(St2, d, u)

        # Current element back to the indeterminate
        St.indeterminate.append(chi)

