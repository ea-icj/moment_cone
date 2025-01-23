from dataclasses import dataclass

import itertools

from .typing import *
from .utils import short_prod
from .group import *
from .tau import *
from .weight import Weight
from .rep import *
from .rings import matrix, Matrix, ZZ

#from typing import TypeVar
#T=TypeVar('T')
#from sage.all import matrix


__all__ = (
    "find_hyperplanes_reg_mod_outer",
    "find_hyperplanes_reg_impl",
    "check_hyperplane_dim",
    "has_too_much_geq_weights",
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
                    S_positive: list[Weight],
                    sym : list[int]=None) -> None:
    """ Determining the sign of weight from S_input by comparing it to the curring weight chi """
    idx = 0
    while idx < len(S_input): # Loop on S_input with while since S_input changes
        chi2 = S_input[idx]
        if chi2.leq(chi,sym) : #chi2<=chi : #if : #if chi<=chi2:# TODO : switch leq
            S_negative.append(chi2)
            smart_remove(S_input, idx)
        elif chi.leq(chi2,sym) : #chi.leq(chi2,sym):
            S_positive.append(chi2)
            smart_remove(S_input, idx)
        else:
            idx += 1 # index incremented only when not element where removed

def check_hyperplane_dim(weights: Sequence[Weight], exp_dim: int) -> bool:
    """ 
    Check that the given sequence of weights span a linear subspace of dimension exp_dim in X^*(T)_QQ
    """
    
    if len(weights) < exp_dim:
        return False

    v=weights[0]
    M = matrix([v.as_vector.list() for v in weights])
    return M.rank(algorithm="flint") == exp_dim # Flint algorithm is faster than the default one

def has_too_much_geq_weights(chi: Weight, weights: list[Weight], V: Representation, u: int, sym: list[int]) -> bool:
    """ 
    True for weights with more that u weights bigger for the order leq 
    """
    if V.type=='kron':
        leq_cnt = short_prod(c + 1 for c in chi.as_list) - 1 #TODO : short est inutile ici
        return leq_cnt > u
    else : 
        leq_cnt=0
        for chi1 in weights:
            if chi.leq(chi1,sym) : # Vérifier que l'ordre est le bon
                leq_cnt+=chi1.mult
        leq_cnt-=chi.mult
        return leq_cnt > u
        

def find_hyperplanes_reg_mod_outer(weights: list[Weight], V: Representation, u: int, sym:list[int]=None) -> Iterable[list[Weight]]:
    """
    Returns the subsets of weights, each set generating an hyperplane in X^*(T) likely to be the orthogonal of a dominant 1-parameter subgroup tau, such that there is at most u weights we of V with tau(we)>0

    Example:
    >>> from cone import *
    >>> d = Dimension((4, 4, 4))
    >>> hp = list(find_hyperplanes_reg(d, 4**3))
    >>> print("Number of raw hyperplanes:", len(hp))
    Number of raw hyperplanes: 3622
    """
    exp_dim=V.dim_cone-1
    St = WeightSieve([], [], [], [], [])
    
    for chi in weights:
        if has_too_much_geq_weights(chi, weights, V, u, sym):
            St.negative.append(chi)
        else:
            St.indeterminate.append(chi)

    if V.type != 'kron': #Trivial outer
        yield from find_hyperplanes_reg_impl(St, V.G, u,exp_dim,sym)
                                                 

    else :
      for chi in V.weights_mod_outer:
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
        yield from find_hyperplanes_reg_impl(St2, V.G, u,exp_dim)

        # Removing symmetries
        for chi2 in chi.orbit_symmetries(V.G.outer):
            St.indeterminate.remove(chi2)
            St.excluded.append(chi2)


def find_hyperplanes_reg_impl(St: WeightSieve, G: LinGroup, u: int,exp_dim: int, sym:list[int]=None) -> Iterable[list[Weight]]:
    """ 
    Recursive part to find the hyperplane candidates.
    u is the maximal number of positive weights
    exp_dim is the expected dimension of the span of the weights.
    sym is a data encoding the partial order on the weights and the definition of dominant tau.

    St is a list of weights of the maximal torus of G.
    sym determine a parabolic subgroup of the Weyl subgroup of G.
    We consider only 1-PS strictly domninant for sym.
    sym also determines a dominance order on the weights
    The weights have a multiplicity.
    """
  
    if check_hyperplane_dim(St.zero, exp_dim):
        # Candidate hyperplane if the dimension is appropriate
        yield St.zero

    elif len(St.zero) + len(St.indeterminate) >= exp_dim and len(St.indeterminate) > 0:
        # Next element to consider
        chi = St.indeterminate.pop()
        St2 = St.copy()
        #St3 = St.copy()
        
        # Two possible actions with this element:

        # 1. We explore the branch where it is excluded from the possible zero elements
        ### TODO : suppress and replace St3 by St
        
        St.excluded.append(chi)
        St3=St.copy()
        yield from find_hyperplanes_reg_impl(St, G, u,exp_dim,sym) # TODO : St n'est pas modifié par cet appel ?
        St.excluded.pop()

        # 2. We explore the branch where it is defined as a zero element (on the hyperplane)
        
        St2.zero.append(chi)

        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(chi, St2.indeterminate, St2.negative, St2.positive,sym)
        sign_assignment(chi, St2.excluded, St2.negative, St2.positive,sym)

        # TODO : Est-il utile de regarder la condition tau régulier ici ??

        # 2.2 Continuing if there are not too much positive elements
        
        if sum([chi.mult for chi in St2.positive]) <= u:
            yield from find_hyperplanes_reg_impl(St2, G, u,exp_dim,sym)

        # Current element back to the indeterminate
        ## TODO : uncomment
        St.indeterminate.append(chi)



