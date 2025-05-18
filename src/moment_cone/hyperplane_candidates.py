from dataclasses import dataclass
#from numba import njit

import itertools
import numpy as np
from numpy.typing import NDArray

from .typing import *
from .utils import short_prod
from .linear_group import *
from .tau import *
from .weight import *
from .representation import *
from .rings import matrix, Matrix, ZZ



__all__ = (
    "find_hyperplanes_reg_mod_outer",
    "find_hyperplanes_reg_impl",
    "check_hyperplane_dim",
    "has_too_much_geq_weights",
)


@dataclass(slots=True)
class WeightSieve:
    """ Sieve used to categorized weights during the search of hyperplane """
    indeterminate: list[int] # Weights whose status is not currently determinate
    excluded: list[int] # Weights excluded from being zero
    zero: list[int] # Weights considered as zero
    positive: list[int] # Positive weights
    negative: list[int] # Negative weights

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


def sign_assignment(id_chi: int,
                    S_input: list[int],
                    S_negative: list[int],
                    S_positive: list[int],
                    MO: NDArray[np.int8]) -> None:
    """ Determining the sign of weight from S_input by comparing it to the curring weight chi """
    idx = 0
    while idx < len(S_input): # Loop on S_input with while since S_input changes
        j = S_input[idx]
        if MO[id_chi,j] == 1 :
            S_positive.append(j)
            smart_remove(S_input, idx)
            #St[j]=1 # chi2 is positive
        elif MO[j,id_chi] == 1 :
            S_negative.append(j)
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
    rank: int = M.rank(algorithm="flint")
    return rank == exp_dim # Flint algorithm is faster than the default one
        
def has_too_much_geq_weights(chi: Weight, weights: Sequence[Weight], V: Representation, u: int, sym: Optional[Sequence[int]] = None) -> bool:
    """ 
    True for weights with more that u weights bigger for the order leq 
    """
    if isinstance(V, KroneckerRepresentation):
        assert(isinstance(chi, WeightAsList))
        leq_cnt = short_prod(c + 1 for c in chi.as_list) - 1 #TODO : short est inutile ici
        return leq_cnt > u
    else : 
        leq_cnt=0
        for chi1 in weights:
            if chi.leq(chi1,sym) : # Vérifier que l'ordre est le bon
                leq_cnt+=chi1.mult
        leq_cnt-=chi.mult
        return leq_cnt > u
    
def find_hyperplanes_reg_mod_outer(weights: Sequence[Weight], V: Representation, u: int, sym: Optional[Sequence[int]] = None) -> Iterable[list[Weight]]:
    """
    Returns the subsets of weights, each set generating an hyperplane in X^*(T) likely to be the orthogonal of a dominant 1-parameter subgroup tau, such that there is at most u weights we of V with tau(we)>0

    Example:
    >>> from moment_cone import *
    >>> G = LinearGroup((4, 4, 4))
    >>> V = KroneckerRepresentation(G)
    >>> hp = list(find_hyperplanes_reg_mod_outer(V.all_weights, V, 4**3))
    >>> print("Number of raw hyperplanes:", len(hp))
    Number of raw hyperplanes: 1604
    """
    exp_dim=V.dim_cone-1
    
    ## TODO : améliorer has_too_much_geq_weights dans le cas Kronecker pour ordonner les poids de V.weights_mod_outer par #{\chi'>\chi} décroissant
    # La fonction qui calcule ce cardinal est leq_cnt = short_prod(c + 1 for c in chi.as_list) - 1

    ## TODO : pour parallelise !='kron' il suffit de faire la même boucle même si modulo outer est trivial
    # On peut même générer les weights directement ordonnés de cette manière
    
    # We cancel weights that are necessarily negative
    weights_free = [chi for chi in weights if not has_too_much_geq_weights(chi, weights, V, u, sym)]
    weights_free_mod_outer = [weights_free.index(chi) for chi in V.weights_mod_outer if chi in weights_free]
    
    # Orbit as a dictionary
    orbit_as_dic_idx: Optional[dict[int, list[int]]] = None
    if isinstance(V, KroneckerRepresentation):
        orbit_as_dic_idx = {i: [weights_free.index(chi2) for chi2 in chi.orbit_symmetries(V.G.outer)] for i, chi in enumerate(weights_free)}

    # Preparatory: compute the incidence matrix of the dominance order
    #               and put the multiplicities in an nparray
    #               we use the index in weights_free
    dom_order_matrix = np.zeros((len(weights_free),len(weights_free)), dtype=np.int8)
    mult_chi_tab = np.zeros((len(weights_free),), dtype=np.int8)
    for i, chi1 in enumerate(weights_free):
        mult_chi_tab[i] = chi1.mult
        for j, chi2 in enumerate(weights_free):
            if chi1.leq(chi2,sym) :
                dom_order_matrix[i,j]=1
            elif chi2.leq(chi1,sym) :
                dom_order_matrix[i,j]=-1   


    # Initialisation of St
    St = WeightSieve([], [], [], [], [])
    St.indeterminate=[i for i in range(len(weights_free))]
    
    if isinstance(V, ParticleRepresentation): #Trivial outer
        raise NotImplementedError
        yield from find_hyperplanes_reg_impl(
            weights_free,
            mult_chi_tab,
            St, V.G, u,exp_dim,dom_order_matrix, orbit_as_dic_idx)
                                                 

    else : # Kronecker
        #for chi in V.weights_mod_outer:
        from itertools import product
        N: int = 2
        for choices in product((False, True), repeat=N):
            yield from find_hyperplanes_reg_inner(
                choices,
                weights_free,
                V,
                u,
                exp_dim,
                mult_chi_tab,
                dom_order_matrix,
                weights_free_mod_outer,
                orbit_as_dic_idx,
            )


def find_hyperplanes_reg_inner(
        choices: tuple[bool, ...], #: chosen branch
        weights_free: Sequence[Weight],
        V: Representation,
        u: int,
        exp_dim: int,
        MW: NDArray[np.int8],
        MO: NDArray[np.int8],
        weights_free_mod_outer: list[int],
        orbit_as_dic_idx: Optional[dict[int, list[int]]] = None,
        ) -> Iterable[list[Weight]]:
    
    St = init_sieve_for_family(choices, weights_free, V, MO, weights_free_mod_outer, orbit_as_dic_idx)
    return find_hyperplanes_reg_impl(weights_free, MW, St, u, exp_dim, MO)


def init_sieve_for_family(choices: tuple[bool, ...], weights_free: Sequence[Weight], V: Representation, MO: NDArray[np.int8], weights_free_mod_outer: list[int], orbit_as_dic_idx: Optional[dict[int, list[int]]] = None) -> WeightSieve:
    St = WeightSieve([], [], [], [], [])
    St.indeterminate=[i for i in range(len(weights_free))]
    if isinstance(V, ParticleRepresentation): 
        for j, test in enumerate(choices):
            if len(St.indeterminate) == 0: break
            # Checking if the element is indeterminate
            id_chi = St.indeterminate.pop()
            if test:
                # We put chi in the hyperplane
                St.zero.append(id_chi)

                # Deducing sign of lower and upper elements
                sign_assignment(id_chi, St.indeterminate, St.negative, St.positive, MO)
                sign_assignment(id_chi, St.excluded, St.negative, St.positive, MO)
            else:
                # We exclude chi from the hyperplane
                St.excluded.append(id_chi)
        return St    
        
    else: # KroneckerRepresentation
        assert isinstance(V, KroneckerRepresentation) and orbit_as_dic_idx is not None
        N = len(choices)
        try:
            first_true = choices.index(True)  
        except ValueError:
            first_true = N

        #i=0
        #for chi in V.weights_mod_outer:
        for id_chi in weights_free_mod_outer[:first_true]: # The False at the beginning    
            #if i>=first_True: break
            #try:
            #    id_chi=weights.index(chi)  # meaning indeterminante 
            #except ValueError:
            #    continue
            #for chi2 in chi.orbit_symmetries(V.G.outer):
            for id_chi2 in orbit_as_dic_idx[id_chi]:        
                St.indeterminate.remove(id_chi2)
                St.excluded.append(id_chi2)
            #i+=1

        if  first_true == N or len(weights_free_mod_outer) <= first_true  : # Only false or no more free weights #next_chi == None :
            return St
        
        next_id_chi = weights_free_mod_outer[first_true]
        #id_chi=weights_free.index(next_chi)        
        # Put next_chi in the hyperplane
        St.zero.append(next_id_chi)# put chi in the hyperplane
        smart_remove(St.indeterminate, next_id_chi)

        # Deducing sign of lower and upper elements
        sign_assignment(next_id_chi, St.indeterminate, St.negative, St.positive, MO)
        sign_assignment(next_id_chi, St.excluded, St.negative, St.positive, MO)
                           
        for j, test in enumerate(choices[first_true+1:]):
            if len(St.indeterminate) == 0: break
            # Checking if the element is indeterminate
            id_chi = St.indeterminate.pop()
            if test:
                # We put chi in the hyperplane
                St.zero.append(id_chi)

                # Deducing sign of lower and upper elements
                sign_assignment(id_chi, St.indeterminate, St.negative, St.positive,MO)
                sign_assignment(id_chi, St.excluded, St.negative, St.positive,MO)
            else:
                # We exclude chi from the hyperplane
                St.excluded.append(id_chi)
        return St    
        
            


def find_hyperplanes_reg_impl(weights: Sequence[Weight],MW: NDArray[np.int8], St: WeightSieve, u: int, exp_dim: int, MO: NDArray[np.int8]) -> Iterable[list[Weight]]:
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
    
    if len(St.zero) >= exp_dim and check_hyperplane_dim([weights[i] for i in St.zero], exp_dim):
        # Candidate hyperplane if the dimension is appropriate
        yield [weights[i] for i in St.zero]

    elif len(St.zero) + len(St.indeterminate) >= exp_dim and len(St.indeterminate) > 0:
        # Next element to consider
        id_chi = St.indeterminate.pop()
        St2 = St.copy()
        
        # Two possible actions with this element:

        # 1. We explore the branch where it is excluded from the possible zero elements
        
        St.excluded.append(id_chi)
        St2.zero.append(id_chi)
        
        yield from find_hyperplanes_reg_impl(weights, MW, St, u, exp_dim, MO)

        # 2. We explore the branch where it is defined as a zero element (on the hyperplane)
        

        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(id_chi, St2.indeterminate, St2.negative, St2.positive,MO)
        sign_assignment(id_chi, St2.excluded, St2.negative, St2.positive,MO)

        # 2.2 Continuing if there are not too much positive elements
        
        if len(St2.positive) <=u or sum([MW[i] for i in St2.positive]) <=u:
            yield from find_hyperplanes_reg_impl(weights, MW, St2, u, exp_dim, MO)


