from utils import *
from permutation import *



def ListWs_Mod_rec(tau : "Tau",pos : int,C_mod : dict[int, int]) -> list["Permutation"]: # List of tuples [w_pos,...,w_len(d)-1] such that U(w)\isom C_mod as tau-module
    d=tau.d
    if pos==len(d)-1:
        return([[w] for w in ListW_Mod(tau,pos,C_mod)])
    Lpos=ListW_subMod(tau,pos,C_mod) # Candidates of w_pos
    res=[]
    for w in Lpos:
        List_Inv=[[pos]+inv for inv in w.inversions()]
        Mw=grading_dictionary(List_Inv, self.dot_root)
        new_C_mod=quotient_C_Mod(C_mod,Mw)
        Lw=ListWs_Mod_rec(tau,pos+1,new_C_mod)
        res+=[[w]+l for l in Lw]
    return(res)

def ListWs_Mod(tau : "Tau") ->  list["Permutation"]:
    Poids_positive=tau.positive_weights
    C_mod={}
    for x in Poids_positive.keys():
        C_mod[x]=len(Poids_positive[x])
    return(ListWs_Mod_rec(tau,0,C_mod))

def Check_Rank_Tpi(ineq : "Inequality", method: "Method") -> bool :
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
    gr = tau.grading_roots(ineq.inversions) # A vérifier
    for x in sorted(gr.keys(),reverse=True): # Run over the possible values of tau.scalar(root) for root inversion of w
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
               uv=action_op_el(root, v, d)
               for row, chi in enumerate(gw[x]): # List of weights such that tau.scalar(chi)=x 
                   M[row,col]=uv[chi.index_ind(d)]
        rank_M = M.change_ring(ring.fraction_field()).rank()
        if rank_M<M.nrows():
               return(False)
    return(True)	       
      
