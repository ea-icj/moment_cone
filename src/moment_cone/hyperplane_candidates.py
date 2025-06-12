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
    nb_positive: list[int] # Number of Positive weights. List of length 1 used to pass the data by address
    nb_pos_checked: list[int]
    #nb_negative: list[int] # Number of Negative weights. List of length 1 used to pass the data by address

    def copy(self) -> "WeightSieve":
        return WeightSieve(
            self.indeterminate.copy(),
            self.excluded.copy(),
            self.zero.copy(),
            self.nb_positive.copy(),
            self.nb_pos_checked.copy(),
            #self.nb_negative.copy(),
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
    l[idx]=l[-1]
    l.pop()                            

def sign_assignment(id_chi: int,
                    S_ex: list[int],
                    S_ind: list[int],
                    nb_positive: list[int],
                    MO: NDArray[np.int8],
                    MW: NDArray[np.int8],
                    V: Representation) -> None:
    """ Determining the sign of weight from S_input by comparing it to the curring weight chi """
    pos_ex = np.flatnonzero(MO[id_chi,:][np.ix_(S_ex,)])
    pos_ind = np.flatnonzero(MO[id_chi,:][np.ix_(S_ind,)])
    if isinstance(V, KroneckerRepresentation):
        nb_positive[0]+=len(pos_ex)+len(pos_ind)
        for idx in pos_ex[::-1]:
            smart_remove(S_ex,idx)
        for idx in pos_ind[::-1]:
            smart_remove(S_ind,idx) 
    else :
        for idx in pos_ex[::-1]:
            nb_positive[0]+=MW[S_ex[idx]]
            smart_remove(S_ex,idx)
        for idx in pos_ind[::-1]:
            nb_positive[0]+=MW[S_ind[idx]]
            smart_remove(S_ind,idx)             

    neg_ind = np.flatnonzero(MO[:,id_chi][np.ix_(S_ind,)])    
    for idx in neg_ind[::-1]:
        smart_remove(S_ind,idx)
    neg_ex = np.flatnonzero(MO[:,id_chi][np.ix_(S_ex,)])     
    for idx in neg_ex[::-1]:
        smart_remove(S_ex,idx)


def put_negative(S_ex: list[int],
                    S_ind: list[int],
                    current_u: int,
                    MO: NDArray[np.int8],
                    MW: NDArray[np.int8]) -> None:
    """ Determining the sign of weight from S_input by comparing it to the curring weight chi """
    # TODO : faire pour particule
    NbSup = MO[np.ix_(S_ex+S_ind,S_ex+S_ind)].sum(axis=1)
    l_ex=len(S_ex) 
    upu=current_u+1
    
    To_be_deleted_ex = np.flatnonzero(NbSup[:l_ex] > upu)
    To_be_deleted_ind = np.flatnonzero(NbSup[l_ex:] > upu)

    for idx in To_be_deleted_ex[::-1]:
        smart_remove(S_ex,idx) 
    for idx in To_be_deleted_ind[::-1]:
        smart_remove(S_ind,idx)     

def check_hyperplane_dim(weights: Sequence[Weight], exp_dim: int) -> bool:
    """ 
    Check that the given sequence of weights span a linear subspace of dimension exp_dim in X^*(T)_QQ
    """  
    if len(weights) < exp_dim:
        return False

    M = matrix([v.as_vector.list() for v in weights])
    rank: int = M.rank(algorithm="flint") # Flint algorithm is faster than the default one
    return rank == exp_dim 
        
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
            if chi.leq(chi1,sym) : 
                leq_cnt+=chi1.mult
        leq_cnt-=chi.mult
        return leq_cnt > u
    
def sort_criterium(id_chi: int, weights_free: list[int],V: Representation,MO: NDArray[np.int8]) -> tuple[int, int]:
    nb_sup = MO[id_chi, :].sum() - 1
    nb_inf = MO[:, id_chi].sum() - 1
    nb_orbit=len(list(weights_free[id_chi].orbit_symmetries(V.G.outer)))
    return nb_orbit, -1.8 * nb_sup - nb_inf

def best_index_for_sign(L: list[int],MO: NDArray[np.int8],coeff: float) -> tuple[int, int]:
    """
    Return the index of an element of L such that nb_inf + nb_sup is maximal
    """
    A = MO[np.ix_(L, L)]
    scores = coeff * A.sum(axis=1) + A.sum(axis=0) # coef 2 to insist on positive weights because of u
    best_i = int(np.argmax(scores))
    return best_i, L[best_i]

def best_index_for_sign2(S_ind: list[int],S_ex: list[int],MO: NDArray[np.int8],coeff: int) -> tuple[int]:
    """
    Return the index of an element of L such that nb_inf + nb_sup is maximal
    """
    A = MO[np.ix_(S_ind, S_ind)]
    if coeff==1 : 
        scores = coeff * A.sum(axis=1) + A.sum(axis=0) # coef 2 to insist on poitive weights because of u
    else : 
        B = MO[np.ix_(S_ind, S_ex)]
        scores = coeff * A.sum(axis=1) + A.sum(axis=0)+ 0.1*B.sum(axis=1)  
    #scores = 2 * MO[L][:, L].sum(axis=1) + MO[L][:, L].sum(axis=0)
    #scores = MO[L][:, L].sum(axis=1) + MO[L][:, L].sum(axis=0)
    best_i = int(np.argmax(scores))
    return best_i,S_ind[best_i]


def find_hyperplanes_reg_mod_outer(
        weights: Sequence[Weight],
        V: Representation,
        umax: int,
        sym: Optional[Sequence[int]] = None,
        flatten_level: int = 1, #: Flatten search graph up to given level
    ) -> Iterable[Tau] : #Iterable[list[Weight]]:
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

    exp_dim = V.dim_cone - 1
    #print(V,umax)
    # We cancel weights that are necessarily negative
    weights_almost_free = [chi for chi in weights if not has_too_much_geq_weights(chi, weights, V, umax, sym)]
    
    weights_free: list[Weight] = []
    for chi1 in weights_almost_free:
        nb_sup=len([chi2 for chi2 in weights_almost_free if chi2.leq(chi1,sym)]) # Todo : ou le contraire
        if nb_sup+exp_dim<=len(weights_almost_free)+1 :
            weights_free.append(chi1)
    u=umax-len(weights_almost_free)+len(weights_free)

    #Preparatory: compute the incidence matrix of the dominance order
    #               and put the multiplicities in an nparray
    #               we use the index in weights_free
    dom_order_matrix = np.zeros((len(weights_free),len(weights_free)), dtype=np.int8)
    mult_chi_tab = np.zeros((len(weights_free),), dtype=np.int8)
    for i, chi1 in enumerate(weights_free):
        mult_chi_tab[i] = chi1.mult
        for j, chi2 in enumerate(weights_free):
            if chi1.leq(chi2,sym) :
                dom_order_matrix[i,j]=1
             
    #n = len(weights)
    #from sage.all import Matrix, GF
    #MO = Matrix(GF(2), n, n, False)
    #for i, chi1 in enumerate(weights):
    #    for j, chi2 in enumerate(weights):
    #        if chi1.leq(chi2,sym)  and i != j:
    #            MO[i,j]=1
    #from sage.all import Poset, DiGraph 
    #Pos = Poset(DiGraph(MO))
    #print('Antichaines',V)
    #print(Pos.width())
    #print(len(Pos.maximal_antichains()))


    #print('V',V,'we have ',len(weights_free),'free weoghts over',len(weights))
    
    weights_free_mod_outer = [weights_free.index(chi) for chi in V.weights_mod_outer if chi in weights_free]
    # For parallel : plus lent mais mieux équilibré
    weights_free_mod_outer.sort(key=lambda id_chi: sort_criterium(id_chi, weights_free,V,dom_order_matrix),reverse=True)
           

    # Orbit as a dictionary
    orbit_as_dic_idx: Optional[dict[int, list[int]]] = None
    if isinstance(V, KroneckerRepresentation):
        orbit_as_dic_idx = {i: [weights_free.index(chi2) for chi2 in weights_free[i].orbit_symmetries(V.G.outer)] for i in weights_free_mod_outer}


    # Sort the weights
    
    if isinstance(V, ParticleRepresentation):
        Choices=[tuple([False]*nb_false+[True]) for nb_false in range(min(flatten_level-1,len(weights_free)))]
    else :    
        Choices=[tuple([False]*nb_false+[True]) for nb_false in range(min(flatten_level-1,len(weights_free_mod_outer)))]
    Nb_False=min(flatten_level-1,len(weights_free_mod_outer))    
    new_flatten_level=flatten_level-len(Choices)
    if new_flatten_level>2 :  #Some binary choices can still be done.
        from math import floor
        from itertools import product
        free_choices = new_flatten_level.bit_length() - 1
        b = floor((len(Choices)+1)/2)
        q, r = divmod(free_choices, b)
        Choices = (
            Choices[:r] + 
            [x + p for x in Choices[r:b] for p in product((False, True), repeat=q)] + 
            Choices[b:]
        )
        Choices = (
            [x + p for x in Choices[:r] for p in product((False, True), repeat=q+1)] + 
            Choices[r:]
        )

    # Initialisation of St
    St = WeightSieve([], [], [], [0],[0])
    St.indeterminate=[i for i in range(len(weights_free))]
    
    # Task corresponding to [False]*Nb_False
    if isinstance(V, ParticleRepresentation): #Trivial outer
        for id_chi in weights_free_mod_outer[:Nb_False]: # One more False. 
            St.indeterminate.remove(id_chi)
            St.excluded.append(id_chi)
        yield from find_hyperplanes_reg_impl(
            weights_free,
            V,
            mult_chi_tab,
            St,u,exp_dim,dom_order_matrix,sym=sym)
                                                 
    else : # Kronecker
        assert orbit_as_dic_idx is not None
        for id_chi in weights_free_mod_outer[:Nb_False]: # One more False.            
            for id_chi2 in orbit_as_dic_idx[id_chi]:        
                St.indeterminate.remove(id_chi2)
                St.excluded.append(id_chi2)
        for i in range(len(weights_free_mod_outer[Nb_False:])):
            if len(St.indeterminate)+len(St.excluded)>u-St.nb_positive[0] :
                coeff_best = 1.8
            else :
                coeff_best = 1             
            # If so, we explore the branch where id_chi is  defined as a zero element (on the hyperplane)
            idx,id_chi = best_index_for_sign(St.indeterminate,dom_order_matrix,coeff_best)
            St2 = St.copy()
            St2.zero.append(id_chi)# put chi in the hyperplane
            smart_remove(St2.indeterminate, idx)
            
            # Deducing sign of lower and upper elements
            sign_assignment(id_chi, St2.excluded, St2.indeterminate, St2.nb_positive,dom_order_matrix,mult_chi_tab,V)

            yield from find_hyperplanes_reg_impl(weights_free,V,mult_chi_tab,St2, u,exp_dim,dom_order_matrix)    
            
            # Now we exclude id_chi and its orbit
            for id_chi2 in [lst for lst in orbit_as_dic_idx.values() if id_chi in lst][0]:        
                St.indeterminate.remove(id_chi2)
                St.excluded.append(id_chi2)
    
    # Autres tâches
    for choices in Choices:  # Liste de tâches    
        yield from find_hyperplanes_reg_inner(
            choices,
            weights_free,
            V,
            u,
            exp_dim,
            mult_chi_tab,
            dom_order_matrix,
            weights_free_mod_outer,
            orbit_as_dic_idx=orbit_as_dic_idx,sym=sym
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
        orbit_as_dic_idx: Optional[dict[int, list[int]]] = None, sym: Optional[Sequence[int]] = None
        ) -> Iterable[Tau]: 
    
    St = init_sieve_for_family(choices, weights_free, V, MO, MW, weights_free_mod_outer, orbit_as_dic_idx)
    return find_hyperplanes_reg_impl(weights_free, V, MW, St, u, exp_dim, MO,sym=sym)


def init_sieve_for_family(choices: tuple[bool, ...], weights_free: Sequence[Weight], V: Representation, MO: NDArray[np.int8], MW: NDArray[np.int8], weights_free_mod_outer: list[int], orbit_as_dic_idx: Optional[dict[int, list[int]]] = None) -> WeightSieve:
    St = WeightSieve([], [], [], [0], [0])
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
                sign_assignment(id_chi, St.excluded, St.indeterminate, St.nb_positive, MO,MW,V)
            else:
                # We exclude chi from the hyperplane
                St.excluded.append(id_chi)    
        
    else: # KroneckerRepresentation
        assert isinstance(V, KroneckerRepresentation) and orbit_as_dic_idx is not None
        N = len(choices)
        try:
            first_true = choices.index(True)  
        except ValueError:
            first_true = N

        for id_chi in weights_free_mod_outer[:first_true]: # The False at the beginning    
            for id_chi2 in orbit_as_dic_idx[id_chi]:        
                St.indeterminate.remove(id_chi2)
                St.excluded.append(id_chi2)
            
        if  first_true == N or len(weights_free_mod_outer) <= first_true  : # Only false or no more free weights #next_chi == None :
            return St
        
        next_id_chi = weights_free_mod_outer[first_true]        
        # Put next_chi in the hyperplane
        St.zero.append(next_id_chi)# put chi in the hyperplane
        St.indeterminate.remove(next_id_chi)
        
        # Deducing sign of lower and upper elements
        sign_assignment(next_id_chi, St.excluded, St.indeterminate, St.nb_positive, MO,MW,V)
           
        for j, test in enumerate(choices[first_true+1:]):
            if len(St.indeterminate) == 0: 
                break
            # Checking if the element is indeterminate
            id_chi = St.indeterminate.pop()
            if test:
                # We put chi in the hyperplane
                St.zero.append(id_chi)

                # Deducing sign of lower and upper elements
                sign_assignment(id_chi, St.excluded, St.indeterminate, St.nb_positive,MO,MW,V)
            else:
                # We exclude chi from the hyperplane
                St.excluded.append(id_chi)
    return St    
        
            


def find_hyperplanes_reg_impl(weights: Sequence[Weight],V: Representation,MW: NDArray[np.int8], St: WeightSieve, u: int, exp_dim: int, MO: NDArray[np.int8], sym: Optional[Sequence[int]] = None) -> Iterable[Tau]: # Tau for V.G
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
        # Candidate hyperplane if the dimension is appropriate. Computation of the dominant equation if there exists.
        taured=Tau.from_zero_weights([weights[i] for i in St.zero], V)
        if isinstance(V, KroneckerRepresentation) :
            taured_test_dom=taured
        else:
            assert sym is not None
            taured_test_dom=Tau.from_flatten(taured.flattened,LinearGroup(sym))
        if taured_test_dom.is_dom_reg : # We keep only dominant regular 1-PS
            yield taured
        elif taured_test_dom.opposite.is_dom_reg :
            yield from taured.opposite
        

    elif len(St.zero) + len(St.indeterminate) >= exp_dim and len(St.indeterminate) > 0:
        if len(St.indeterminate)+len(St.excluded)>u-St.nb_positive[0] :
            coeff_best = 1.8
        else :
            coeff_best = 1    
        # Next element to consider
        idx,id_chi = best_index_for_sign(St.indeterminate,MO,coeff_best)
        smart_remove(St.indeterminate, idx)
        St2 = St.copy()
        
        # Two possible actions with this element:

        # 1. We explore the branch where it is excluded from the possible zero elements
        
        St.excluded.append(id_chi)
        St2.zero.append(id_chi)
        
        yield from find_hyperplanes_reg_impl(weights, V, MW, St, u,exp_dim, MO,sym=sym)

        # 2. We explore the branch where it is defined as a zero element (on the hyperplane)
        
        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(id_chi, St2.excluded, St2.indeterminate, St2.nb_positive,MO,MW,V)
        if u>=St2.nb_positive[0]  and St2.nb_pos_checked[0] < St2.nb_positive[0] and  len(St2.indeterminate)+len(St2.excluded)>u-St2.nb_positive[0] : # Otherwise it is unuseful
            St2.nb_pos_checked[0] = St2.nb_positive[0]
            put_negative(St2.excluded, St2.indeterminate, u-St2.nb_positive[0],MO,MW)

        # 2.2 Continuing if there are not too much positive elements
        
        if St2.nb_positive[0] <=u:    
            yield from find_hyperplanes_reg_impl(weights, V, MW, St2, u, exp_dim, MO,sym=sym)    


