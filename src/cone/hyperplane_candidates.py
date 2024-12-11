from dataclasses import dataclass
import itertools

from .typing import *
from .weight import Weight
from .permutation import Permutation
from .dimension import Dimension
from .utils import short_prod
from .tau import *
#from sage.all import matrix, ZZ, QQ # type: ignore
from .rings import matrix, Matrix, ZZ

__all__ = (
    "find_hyperplanes_reg",
    "hyperplane_matrix",
    "find_hyperplanes_reg_mod_sym_dim",
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


def smart_remove(l: list[T], idx: int) -> None:
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
                    S_positive: list[Weight]) -> None:
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

def hyperplane_matrix(S: Sequence[Weight], d: Dimension) -> Matrix:
    """ S is a set of weights. 
    The coordinates of each weight in the basis of fundamental weights are in {0,1}. 
    Each column of the returned matrix is made of the coordinates of a weight in S
    (1st coordinate 1 for C-component, len(d) other blocks of respective size di)
    
    Example :
    >>> d = Dimension((3, 3, 3))
    >>> SW = list(Weight.all(d))[0:5]
    >>> for chi in SW:
    ...     print(chi)
    Weight((0, 0, 0), idx: 0)
    Weight((0, 0, 1), idx: 1)
    Weight((0, 0, 2), idx: 2)
    Weight((0, 1, 0), idx: 3)
    Weight((0, 1, 1), idx: 4)
    >>> hyperplane_matrix(SW, d)
    [1 1 1 1 1]
    [1 1 1 1 1]
    [0 0 0 0 0]
    [0 0 0 0 0]
    [1 1 1 0 0]
    [0 0 0 1 1]
    [0 0 0 0 0]
    [1 0 0 1 0]
    [0 1 0 0 1]
    [0 0 1 0 0]
    """
    M = matrix(ZZ, d.sum + 1, len(S))
    # We could use Blocks but matrix doesn't support len
    shift = tuple(itertools.accumulate(d, initial=1))[:-1]
    for j, chi in enumerate(S):
        M[0, j] = 1
        for c, s in zip(chi, shift):
            M[s + c, j] = 1
    return M

def check_hyperplane_dim(S: Sequence[Weight], d: Dimension) -> bool:
    """ Check that the given sequence of weight has the dimension of an hyperplane in X^*(T)"""
    target_rank = d.sum - len(d)  #Since T lies in C^*\times product of the SL_{d_i}, it has codimension len(d) in \hat T (torus of C^*\times product of the GL_{d_i}
    if len(S) < target_rank:
        return False
    
    M = hyperplane_matrix(S, d)
    return M.rank(algorithm="flint") == target_rank # Flint algorithm is faster than the default one

def has_too_much_geq_weights(chi: Weight, u: int) -> bool:
    """ True for weights with more that u weights bigger for the order leq """
    leq_cnt = short_prod(c + 1 for c in chi) - 1
    # FIXME: verify this function, it differ from the original version
    return leq_cnt > u

def find_hyperplanes_reg(d: Dimension, u: int) -> Iterable[list[Weight]]:
    """
    Returns sets of weights, each set generating an hyperplane in X^*(T) likely to be the orthogonal of a dominant 1-parameter subgroup tau, such that there is at most u weights we of V with tau(we)>0

    Example:
    >>> from cone import *
    >>> d = Dimension((4, 4, 4))
    >>> hp = list(find_hyperplanes_reg(d, 4**3))
    >>> print("Number of raw hyperplanes:", len(hp))
    Number of raw hyperplanes: 3622
    """
    St = WeightSieve([], [], [], [], [])
    for chi in Weight.all(d):
        if has_too_much_geq_weights(chi, u):
            St.negative.append(chi)
        else:
            St.indeterminate.append(chi)
    return find_hyperplanes_reg_impl(St, d, u)

def find_hyperplanes_reg_mod_sym_dim(d: Dimension, u: int) -> Iterable[list[Weight]]:
    """
    Same as find_hyperplanes, with results up to the action of the symmetries of d.
    
    It should be faster since less possibilities are explored.

    Example:
    >>> from cone import *
    >>> d = Dimension((4, 4, 4))
    >>> hp = list(find_hyperplanes_reg_mod_sym_dim(d, 4**3))
    >>> print("Number of raw hyperplanes:", len(hp))
    Number of raw hyperplanes: 1604
    """
    # Earlier filtering
    St = WeightSieve([], [], [], [], [])
    for chi in Weight.all(d):
        if has_too_much_geq_weights(chi, u):
            St.negative.append(chi)
        else:
            St.indeterminate.append(chi)

    # Calling the recursive function modulo the symmetries of d
    for chi in Weight.all_mod_sym_dim(d):
        # Checking if the element is indeterminate
        try:
            idx = St.indeterminate.index(chi)
        except ValueError:
            continue

        # If so, we explore the branch where it is defined as a zero element (on the hyperplane)
        St2 = St.copy()
        St2.zero.append(chi)
        smart_remove(St2.indeterminate, idx)

        # Deducing sign of lower and upper elements
        sign_assignment(chi, St2.indeterminate, St2.negative, St2.positive)
        sign_assignment(chi, St2.excluded, St2.negative, St2.positive)

        # Further exploring the branch
        yield from find_hyperplanes_reg_impl(St2, d, u)

        # Removing symmetries
        for orbit_chi in chi.orbit_symmetries(d.symmetries):
            St.indeterminate.remove(orbit_chi)
            St.excluded.append(orbit_chi)

def find_hyperplanes_reg_impl(St: WeightSieve, d: Dimension, u: int) -> Iterable[list[Weight]]:
    """ Recursive part to find the hyperplane candidates """
    if check_hyperplane_dim(St.zero, d):
        # Candidate hyperplane if the dimension is appropriate
        yield St.zero

    elif len(St.zero) + len(St.indeterminate) >= d.sum - len(d) and len(St.indeterminate) > 0:
        # Next element to consider
        chi = St.indeterminate.pop()

        # Two possible actions with this element:

        # 1. We explore the branch where it is excluded from the possible zero elements
        St.excluded.append(chi)
        yield from find_hyperplanes_reg_impl(St, d, u)
        St.excluded.pop()

        # 2. We explore the branch where it is defined as a zero element (on the hyperplane)
        St2 = St.copy()
        St2.zero.append(chi)

        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(chi, St2.indeterminate, St2.negative, St2.positive)
        sign_assignment(chi, St2.excluded, St2.negative, St2.positive)

        # 2.2 Continuing if there are not too much positive elements
        if len(St2.positive) <= u:
            yield from find_hyperplanes_reg_impl(St2, d, u)

        # Current element back to the indeterminate
        St.indeterminate.append(chi)



