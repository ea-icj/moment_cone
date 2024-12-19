from .typing import *
from .utils import *
from .permutation import *
from .root import *
from .vector_chooser import *
from .rings import matrix
from .spaces import action_op_el
from .tau import Tau
from .inequality import Inequality

__all__ = (
    "ListW_Mod",
    "ListWs_Mod",
    "Check_Rank_Tpi",
)

def ListW_subMod(tau : Tau, pos : int,C_mod : dict[int, int]) -> list[Permutation]:
    "List of permutations w in W^{P(tau[pos])} such that tau.Scalar(Inv(w) in position pos) is a submodule of a C^*-module whose dimension of eigenspaces is encoded by C_mod."
    D=sum(C_mod.values())
    e=tau.d[pos]
    ap = AllPermutationsByLength(e)
    res=[]
    for l in range(min(D+1,int(e*(e-1)/2)+1)):
      for w in ap[l] :
         print('w,tau_mult,minrep',w,tau.reduced.mult[pos],w.is_min_rep(tau.reduced.mult[pos]))
         if w.is_min_rep(tau.reduced.mult[pos]): 
            List_Inv=[Root(pos, *inv) for inv in w.inversions]
            gr=grading_dictionary(List_Inv, tau.dot_root)
            Mw=dictionary_list_lengths(gr)
            if Is_Sub_Mod(Mw,C_mod):
                res.append(w)
    return(res)        

    
def ListW_Mod(tau: Tau,pos : int,C_mod : dict[int, int]) -> list[Permutation]:
    "List of permutations w in W^{P(tau[pos])} such that tau.Scalar(Inv(w) in position pos) is isomorphic to C_mod."
    from .utils import grading_dictionary
    D=sum(C_mod.values())
    e=tau.d[pos]
    ap = AllPermutationsByLength(e)
    res: list[Permutation] = []
    if D>=int(e*(e-1)/2+1):
        return(res)
    for w in ap[D]:
        if w.is_min_rep(tau.reduced.mult[pos]):
            List_Inv=[Root(pos, *inv) for inv in w.inversions]
            gr = grading_dictionary(List_Inv, tau.dot_root)
            Mw=dictionary_list_lengths(gr)
            if Are_Isom_Mod(Mw,C_mod):
                res.append(w)
    return res



def ListWs_Mod_rec(tau: Tau, pos : int, C_mod : dict[int, int]) -> list[list[Permutation]]: # List of tuples [w_pos,...,w_len(d)-1] such that U(w)\isom C_mod as tau-module and w_i\in W^P
    d=tau.d
    print(C_mod)
    if pos==len(d)-1:
        return([[w] for w in ListW_Mod(tau,pos,C_mod)])
    Lpos=ListW_subMod(tau,pos,C_mod) # Candidates of w_pos
    print('pos,Lpos',pos,Lpos)
    res=[]
    for w in Lpos:
        List_Inv=[Root(pos, *inv) for inv in w.inversions]
        gr=grading_dictionary(List_Inv, tau.dot_root)
        Mw=dictionary_list_lengths(gr)
        new_C_mod=quotient_C_Mod(C_mod,Mw)
        Lw=ListWs_Mod_rec(tau,pos+1,new_C_mod)
        res+=[[w]+l for l in Lw]
    print(res)    
    return res

def ListWs_Mod(tau : "Tau") ->  list[list[Permutation]]:
    Poids_positive=tau.positive_weights
    C_mod={}
    for x in Poids_positive.keys():
        C_mod[x]=len(Poids_positive[x])
    print('C_mod départ:',C_mod)    
    return ListWs_Mod_rec(tau,0,C_mod)


def Check_Rank_Tpi(ineq : Inequality, method: Method) -> bool :
    tau=ineq.tau
    d = tau.d # FIXME: Here, d is recreated from scratch, without rings. Should we ensure the uniqueness of the instance of d?

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = d.QI
    elif method == "symbolic":
        ring = d.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    zero_weights = tau.orthogonal_weights
    v = point_vect(zero_weights, d, ring, bounds=(-100, 100))
    gw = tau.grading_weights
    gr = tau.grading_roots_in(ineq.inversions) # A vérifier
    for x in sorted(gr.keys(),reverse=True): # Run over the possible values of tau.scalar(root) for root inversion of w
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
               uv=action_op_el(root, v, d)
               for row, chi in enumerate(gw[x]): # List of weights such that tau.scalar(chi)=x 
                   M[row,col]=uv[chi.index_in(d)]
        rank_M = M.change_ring(ring.fraction_field()).rank()
        if rank_M<M.nrows():
               return(False)
    return(True)	       
      

