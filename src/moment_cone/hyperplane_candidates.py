from dataclasses import dataclass
from flint import fmpz_mat # type: ignore

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
    
    def copy(self) -> "WeightSieve":
        return WeightSieve(
            self.indeterminate.copy(),
            self.excluded.copy(),
            self.zero.copy(),
            self.nb_positive.copy(),
            self.nb_pos_checked.copy(),
        )


def is_parallelizable(St: WeightSieve,nbTrue,nbFalse) -> bool :
    return(len(St.zero)>=nbTrue or len(St.excluded)>=nbFalse)

def smart_remove(l: list[T], idx: int) -> None:
    """
    Remove an element of a list by swapping it with the last element
    idx is the index in l of the element to be removed

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

def check_hyperplane_dim(weights: Sequence[int], exp_dim: int, M_weights : NDArray[np.int8]) -> bool:
    """ 
    Check that the given sequence of weights span a linear subspace of dimension exp_dim in X^*(T)_QQ
    """  
    if len(weights) < exp_dim:
        return False

    M_np = M_weights[:, weights]
    M_flint = fmpz_mat(M_np.tolist())
    return cast(int, M_flint.rank()) == exp_dim  
        
def has_too_much_geq_weights(chi: Weight, weights: Sequence[Weight], V: Representation, u: int, sym: Optional[Sequence[int]] = None) -> bool:
    """ 
    True for weights with more that u weights bigger for the order leq 
    """
    if isinstance(V, KroneckerRepresentation):
        assert(isinstance(chi, WeightAsList))
        leq_cnt = short_prod(c + 1 for c in chi.as_list) - 1 
        return leq_cnt > u
    else : 
        leq_cnt=0
        for chi1 in weights:
            if chi.leq(chi1,sym) : 
                leq_cnt+=chi1.mult
        leq_cnt-=chi.mult
        return leq_cnt > u
    
    
def best_index_for_sign(L: list[int],MO: NDArray[np.int8],coeff: float) -> tuple[int, int]:
    """
    Return the tuple (index,value) of an element of L such that nb_inf + nb_sup is maximal
    """
    A = MO[L][:,L]
    scores = coeff * A.sum(axis=1) + A.sum(axis=0) # coeff to insist on positive weights because of u
    best_i = int(np.argmax(scores))
    return best_i, L[best_i]


def find_hyperplanes_reg_mod_outer(
        weights: Sequence[Weight],
        V: Representation,
        umax: int,
        sym: Optional[Sequence[int]] = None
    ) -> Iterable[Tau] : #Iterable[list[Weight]]:
    """
    Returns the subsets of weights, each set generating an hyperplane in X^*(T) likely to be the orthogonal of a dominant 1-parameter subgroup tau, such that there is at most u weights we of V with tau(we)>0

    Example:
    >>> from moment_cone import *
    >>> G = LinearGroup((4, 4, 4))
    >>> V = KroneckerRepresentation(G)
    >>> hp = list(find_hyperplanes_reg_mod_outer(V.all_weights, V, 4**3))
    >>> print("Number of raw hyperplanes:", len(hp))
    Number of raw hyperplanes: 796
    """

    exp_dim = V.dim_cone - 1

    # We cancel weights that are necessarily negative
    weights_almost_free = [chi for chi in weights if not has_too_much_geq_weights(chi, weights, V, umax, sym)]
    
    # We cancel weights that are necessarily positive
    # and we adjust umax in u accordingly

    weights_free: list[Weight] = []
    for chi1 in weights_almost_free:
        nb_sup=len([chi2 for chi2 in weights_almost_free if chi2.leq(chi1,sym)]) 
        if nb_sup+exp_dim<=len(weights_almost_free)+1 :
            weights_free.append(chi1)
    u=umax-len(weights_almost_free)+len(weights_free)

    

    # Initialisation of the criterium to be parallelizable
    NbTrue = 2
    NbFalse = len(weights_free) //10
    
    #Preparatory: Matrix of weights_free

    M_weights = np.array([v.as_vector.list() for v in weights_free], dtype=np.int8).T 
    
    if isinstance(V, KroneckerRepresentation): # Add End0 conditions
        zero_cols = np.zeros((V.G.rank, len(V.G)-1), dtype=np.int8)

        #M_weights = M_weights.augment(matrix(ZZ, [(len(V.G)-1) * [0] for i in range(V.G.rank)]))
        shift_i = 0
        #shift_j = 0len(weights_free)
        for j, d in enumerate(V.G[:-1]):           
            zero_cols[shift_i + d-1,j] = 1
            shift_i += d
        M_weights = np.hstack([M_weights, zero_cols])

    #Preparatory: compute the incidence matrix of the dominance order
    #               and put the multiplicities in an nparray
    #               we use the indeces in weights_free

    dom_order_matrix = np.zeros((len(weights_free),len(weights_free)), dtype=np.int8)
    mult_chi_tab = np.zeros((len(weights_free),), dtype=np.int8)
    for i, chi1 in enumerate(weights_free):
        mult_chi_tab[i] = chi1.mult
        for j, chi2 in enumerate(weights_free):
            if chi1.leq(chi2,sym) :
                dom_order_matrix[i,j]=1

    # List of weights (their indices in weights_free) modulo outer   
    weights_free_mod_outer = [weights_free.index(chi) for chi in V.weights_mod_outer if chi in weights_free]

    # Orbit as a dictionary using indices in weights_free
    orbit_as_dic_idx: Optional[dict[int, list[int]]] = None
    if isinstance(V, KroneckerRepresentation):
        List_orbits=[[weights_free.index(chi2) for chi2 in weights_free[i].orbit_symmetries(V.G.outer)] for i in weights_free_mod_outer]
        orbit_as_dic_idx = {i: orbit for orbit in List_orbits for i in orbit}
    else :
        orbit_as_dic_idx = {i: [i] for i in weights_free_mod_outer}    


    # Initialisation of St
    St = WeightSieve([], [], [], [0],[0])
    # Some weights that are necessarily excluded: if there are not enough incomparable weights to span an hyperplane.
    Excluded = []
    Indeterminate=[]
    for id_chi1, chi1 in enumerate(weights_free):
        nb_comp=len([chi2 for chi2 in weights_free if chi2.leq(chi1,sym) or chi1.leq(chi2,sym)])
        if len(weights_free)-nb_comp<exp_dim-1 :
            St.excluded.append(id_chi1)
        else :
            St.indeterminate.append(id_chi1)    

    #St.indeterminate=[i for i in range(len(weights_free))]

    
    yield from find_hyperplanes_reg_impl(weights_free,V,mult_chi_tab,St, u,exp_dim,dom_order_matrix, M_weights,orbit_as_dic_idx, NbTrue,NbFalse,sym=sym) 



def find_hyperplanes_reg_impl(weights: Sequence[Weight],V: Representation,MW: NDArray[np.int8], St: WeightSieve, u: int, exp_dim: int, MO: NDArray[np.int8], M_weights : NDArray[np.int8], orbit_as_dic_idx:dict[int, list[int]], NbTrue: int,NbFalse: int, sym: Optional[Sequence[int]] = None) -> Iterable[Tau]: # Tau for V.G
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
    
    # Case when St.zero determines an hyperplane
    if len(St.zero) >= exp_dim and check_hyperplane_dim(St.zero, exp_dim, M_weights):
        # Candidate hyperplane if the dimension is appropriate. Computation of the dominant equation if there exists.
        taured=Tau.from_zero_weights(St.zero, M_weights, V)

        if isinstance(V, KroneckerRepresentation) :
            taured_test_dom=taured
        else:
            assert sym is not None
            taured_test_dom=Tau.from_flatten(taured.flattened,LinearGroup(sym))
        if taured_test_dom.is_dom_reg : # We keep only dominant regular 1-PS
            yield taured
        elif taured_test_dom.opposite.is_dom_reg :
            yield taured.opposite
        
    # Case when St.indeterminate is sufficiently big to get hyperplanes
    elif len(St.zero) + len(St.indeterminate) >= exp_dim and len(St.indeterminate) > 0:
        # check if u is still a restrictive condition
        if len(St.indeterminate)+len(St.excluded)>u-St.nb_positive[0] :
            coeff_best = 1.8
        else :
            coeff_best = 1  

        # Next element to consider
        idx,id_chi = best_index_for_sign(St.indeterminate,MO,coeff_best)
        
        St2 = St.copy()
        
        # Two possible actions with this element:

        # 1. We explore the branch where id_chi is excluded from the possible zero elements
        # If no weight is still selected we exclude the orbit of id_chi1
        if len(St.zero) > 0 :
            smart_remove(St.indeterminate, idx)
            St.excluded.append(id_chi)
        else : 
            for id_chi2 in orbit_as_dic_idx[id_chi]:#coucou
                #smart_remove(St.indeterminate, id_chi2)        
                St.indeterminate.remove(id_chi2)
                St.excluded.append(id_chi2)

        if is_parallelizable(St,NbTrue,NbFalse) : #parallel
            yield from find_hyperplanes_reg_impl(weights, V, MW, St, u,exp_dim, MO, M_weights,orbit_as_dic_idx,NbTrue,NbFalse, sym=sym)
        else : # sequential 
            yield from find_hyperplanes_reg_impl(weights, V, MW, St, u,exp_dim, MO, M_weights,orbit_as_dic_idx,NbTrue,NbFalse, sym=sym)

        # 2. We explore the branch where it is defined as a zero element (on the hyperplane)
        St2.zero.append(id_chi)
        smart_remove(St2.indeterminate, idx)
        
        # 2.1 Deducing sign of lower and upper elements
        sign_assignment(id_chi, St2.excluded, St2.indeterminate, St2.nb_positive,MO,MW,V)
        if u>=St2.nb_positive[0]  and St2.nb_pos_checked[0] < St2.nb_positive[0] and  len(St2.indeterminate)+len(St2.excluded)>u-St2.nb_positive[0] : # Otherwise it is unuseful
            St2.nb_pos_checked[0] = St2.nb_positive[0]
            put_negative(St2.excluded, St2.indeterminate, u-St2.nb_positive[0],MO,MW)

        # 2.2 Continuing if there are not too much positive elements
        
        if St2.nb_positive[0] <=u:
            if is_parallelizable(St,NbTrue,NbFalse) : #parallel
                yield from find_hyperplanes_reg_impl(weights, V, MW, St2, u,exp_dim, MO, M_weights,orbit_as_dic_idx,NbTrue,NbFalse, sym=sym)
            else : # sequential 
                yield from find_hyperplanes_reg_impl(weights, V, MW, St2, u,exp_dim, MO, M_weights,orbit_as_dic_idx,NbTrue,NbFalse, sym=sym)    
            


