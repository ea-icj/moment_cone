from dataclasses import dataclass

import itertools
import numpy as np

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

def sign_assignment(id_chi: int,
                    St,
                    MO) -> None:
    """ Determining the sign of weight from S_input by comparing it to the curring weight chi """

    #id_chi = weights.index(chi)
    for j in range(St.size):
        if St[j]>1 : # chi2 is excluded or indeterminate
            if MO[id_chi,j] == 1 : St[j]=1 # chi2 is positive
            elif MO[j,id_chi] == 1 : St[j]=-1 # chi2 is negative
    

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

def has_too_much_geq_weights(chi: Weight, weights: Sequence[Weight], V: Representation, u: int, MO = None) -> bool:
    """ 
    True for weights with more that u weights bigger for the order leq 
    """
    if isinstance(V, KroneckerRepresentation):
        assert(isinstance(chi, WeightAsList))
        leq_cnt = short_prod(c + 1 for c in chi.as_list) - 1 #TODO : short est inutile ici
        return leq_cnt > u
    else : 
        leq_cnt=0
        i=weights.index(chi)
        for j,chi1 in enumerate(weights):
            if MO[i,j]==1 : #chi.leq(chi1,sym) :
                leq_cnt+=chi1.mult
        leq_cnt-=chi.mult
        return leq_cnt > u
        

def find_hyperplanes_reg_mod_outer(weights: Sequence[Weight], V: Representation, u: int, sym: Optional[Sequence[int]] = None) -> Iterable[list[Weight]]:
    """
    Returns the subsets of weights, each set generating an hyperplane in X^*(T) likely to be the orthogonal of a dominant 1-parameter subgroup tau, such that there is at most u weights we of V with tau(we)>0

    Example:
    >>> from cone import *
    >>> G = LinearGroup((4, 4, 4))
    >>> V = KroneckerRepresentation(G)
    >>> hp = list(find_hyperplanes_reg_mod_outer(V.all_weights, V, 4**3))
    >>> print("Number of raw hyperplanes:", len(hp))
    Number of raw hyperplanes: 1604
    """
    exp_dim=V.dim_cone-1
    #St = WeightSieve([], [], [], [], [])
    
    ## TODO : améliorer has_too_much_geq_weights dans le cas Kronecker pour ordonner les poids de V.weights_mod_outer par #{\chi'>\chi} décroissant
    # La fonction qui calcule ce cardinal est leq_cnt = short_prod(c + 1 for c in chi.as_list) - 1

    ## TODO : pour parallelise !='kron' il suffit de faire la même boucle même si modulo outer est trivial
    # On peut même générer les weights directement ordonnés de cette manière
    
    # Preperatory 1 : compute the incidence matrix of the dominance order
    #               and put the multiplicities in an nparray
    dom_order_matrix = np.zeros((len(weights),len(weights)), dtype=np.int8)
    mult_chi_tab = np.zeros((len(weights),), dtype=np.int8)
    for i,chi1 in enumerate(weights):
        mult_chi_tab[i] = chi1.mult
        for j,chi2 in enumerate(weights):
            if chi1.leq(chi2,sym) :
                dom_order_matrix[i,j]=1
            elif chi2.leq(chi1,sym) :
                dom_order_matrix[i,j]=-1   

    St = np.full(len(weights), 7, dtype=np.int8)
    for i,chi in enumerate(weights):
        if has_too_much_geq_weights(chi, weights, V, u, dom_order_matrix):
            St[i]=-1 # chi has to be  negative
        
    if isinstance(V, ParticleRepresentation): #Trivial outer
        yield from find_hyperplanes_reg_impl(weights,mult_chi_tab,St, V.G, u,exp_dim,dom_order_matrix)
                                                 

    else : # Kronecker
      for chi in V.weights_mod_outer: 
        # Checking if the element is indeterminate
        
        id_chi=weights.index(chi)
        try:
            St[id_chi] == 7  # meaning indeterminante 
        except ValueError:
            continue

        # If so, we explore the branch where it is defined as a zero element (on the hyperplane)
        St2 = St.copy()
        St2[id_chi] = 0  # put chi in the hyperplane
        #smart_remove(St2.indeterminate, id_chi)

        # Deducing sign of lower and upper elements
        #sign_assignment(chi, St2.indeterminate, St2.negative, St2.positive)
        sign_assignment(id_chi, St2,dom_order_matrix)
        #sign_assignment(chi, St2.excluded, St2.negative, St2.positive)

        # Further exploring the branch
        yield from find_hyperplanes_reg_impl(weights,mult_chi_tab,St2, V.G, u,exp_dim,dom_order_matrix)

        # Removing symmetries
        for chi2 in chi.orbit_symmetries(V.G.outer):
            St[weights.index(chi2)]=2
            #St.indeterminate.remove(chi2)
            #St.excluded.append(chi2)


def find_hyperplanes_reg_impl(weights: Sequence[Weight],MW, St, G: LinearGroup, u: int, exp_dim: int, MO) -> Iterable[list[Weight]]:
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
    
    
    Nb_indet = np.sum(St == 7)
    Nb_zero = np.sum(St == 0)
    if Nb_zero >= exp_dim :
        List_weights_hyp=[chi for i,chi in enumerate(weights) if St[i]==0]
        if check_hyperplane_dim(List_weights_hyp, exp_dim):
            # Candidate hyperplane if the dimension is appropriate
            yield List_weights_hyp

    elif Nb_zero + Nb_indet >= exp_dim and Nb_indet > 0:
        # Next element to consider
        id_chi = np.argmax(St == 7) # First indeterminate weight
        chi = weights[id_chi]
        #chi = St.indeterminate.pop()
        St2 = St.copy()
        
        # Two possible actions with this element:

        # 1. We explore the branch where it is excluded from the possible zero elements
        
        #St.excluded.append(chi)
        St[id_chi]=2 # exclude chi
        St2[id_chi]=0 # put chi in the hyperplane
        yield from find_hyperplanes_reg_impl(weights,MW,St, G, u,exp_dim,MO)
        #St.excluded.pop()

        # 2. We explore the branch where it is defined as a zero element (on the hyperplane)
        
        #St2.zero.append(chi)

        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(id_chi,St2,MO)
        #sign_assignment(chi, St2,MO)

        # TODO : Est-il utile de regarder la condition tau régulier ici ??

        # 2.2 Continuing if there are not too much positive elements
        
        #if sum([chi.mult for i,chi in enumerate(weights) if St2[i]==1]) <= u:
        if MW[St == 1].sum()<=u:
            yield from find_hyperplanes_reg_impl(weights,MW,St2, G, u,exp_dim,MO)

        # Current element back to the indeterminate
        St[id_chi]=7 
        #St.indeterminate.append(chi)



