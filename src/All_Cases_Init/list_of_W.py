
from .typing import *
from .utils import *
from .permutation import *
from .root import *
from .vector_chooser import *
from .rings import matrix
from .tau import Tau
from .rep import *
from .inequality import Inequality
import operator

__all__ = (
    "ListW_Mod",
    "ListWs_Mod",
    "Check_Rank_Tpi",
)

def ListW_Mod(tau : Tau, pos : int, C_mod : dict[int, int],relation) -> list[Permutation]:
    """
    List of permutations w in W^{P(tau[pos])} such that tau.Scalar(Inv(w) in position pos) satisfies relation (leq or eq) 
    with the C^*-module whose dimension of eigenspaces is encoded by C_mod.
    """
    D=sum(C_mod.values()) # Dimension of the C^*-module
    e=tau.G[pos] # Rank of the current GL
    res=[]
    # Lenghts to run over
    if relation==operator.eq and 2*D> e*(e-1):
        return []
    for w in Permutation.all_min_rep(tau.reduced.mult[pos]):
        if relation(w.length,D): # The length makes possible the expected relation
            List_Inv=[Root(pos, *inv) for inv in w.inversions]
            gr=grading_dictionary(List_Inv, tau.dot_root)
            Mw=dictionary_list_lengths(gr)
            if compare_C_Mod(Mw,C_mod,relation):
                res.append(w)
    return res


def ListWs_Mod_rec(tau: Tau, pos : int, C_mod : dict[int, int]) -> list[list[Permutation]]:
    """ 
    List of tuples [w_pos,...,w_len(d)-1] such that U(w)\isom C_mod as tau-module and w_i\in W^P
    """
    
    G=tau.G
    if pos==len(G)-1: #Only one w has to be find
        relation = operator.eq
        return([[w] for w in ListW_Mod(tau,pos,C_mod,relation)])
    
    Lpos=ListW_Mod(tau,pos,C_mod,operator.le) # Candidates of w_pos
    res=[]
    for w in Lpos:
        List_Inv=[Root(pos, *inv) for inv in w.inversions]
        gr=grading_dictionary(List_Inv, tau.dot_root)
        Mw=dictionary_list_lengths(gr)
        new_C_mod=quotient_C_Mod(C_mod,Mw)
        Lw=ListWs_Mod_rec(tau,pos+1,new_C_mod)
        res+=[[w]+l for l in Lw]
    return res

def ListWs_Mod(tau : Tau,V: Representation) ->  list[list[Permutation]]:
    """
    Initialisation and use the recursive function
    """
    Poids_positive=tau.positive_weights(V)
    C_mod: dict[int, int] = {}
    for x in Poids_positive.keys():
        C_mod[x]=len(Poids_positive[x])
    return ListWs_Mod_rec(tau,0,C_mod)


def Check_Rank_Tpi(ineq : Inequality, V: Representation, method: Method) -> bool :
    """
    Check if Tpi is inversible at a general point of V^\tau. 
    General means randon if the method is probabilist and formal if the method is symbolic.
    The matrix being block trinagular, the function check successively the diagonal blocks.
    """
    tau=ineq.tau
    G = tau.G # FIXME: Here, d is recreated from scratch, without rings. Should we ensure the uniqueness of the instance of d?

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = V.QI
    elif method == "symbolic":
        ring = V.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    zero_weights = tau.orthogonal_weights(V)
    v = point_vect(zero_weights, V, ring, bounds=(-100, 100))
    gw = tau.grading_weights(V)
    gr = tau.grading_roots_in(ineq.inversions) # A vérifier
    for x in sorted(gr.keys(),reverse=True): # Run over the possible values of tau.scalar(root) for root inversion of w
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
               uv=V.action_op_el(root,v)
               for row, chi in enumerate(gw[x]): # List of weights such that tau.scalar(chi)=x 
                   M[row,col]=uv[chi.idx(V)]
        rank_M = M.change_ring(ring.fraction_field()).rank()
        if rank_M<M.nrows():
               return(False)
    return(True)	       


