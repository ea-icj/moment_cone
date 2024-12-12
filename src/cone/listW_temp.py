

def quotient_C_Mod(M1 : dict[int, int], M2 : dict[int, int]) -> dict[int]:
    M = {}
    for p in M1.keys():
        # Valeur par défaut pour M2[p] : 0
        difference = M1[p] - M2.get(p, 0)
        # Ajouter à M uniquement si la différence est non nulle
        if difference != 0:
            M[p] = difference
    return M
    
def ListW_Mod(tau : "Tau",pos : int,C_mod : dict[int, int]) -> list["Permutation"]:
    "List of permutations w such that tau.Scalar(Inv(w) in position pos) is isomorphic to C_mod."
    D=sum(C_mod.values())
    e=tau.d[pos]
    ap = AllPermutationsByLength(e)
    res=[]
    for w in ap[D]:
        List_Inv=[[pos]+inv for inv in w.inversions]
        Mw=grading_dictionary(List_Inv, self.dot_root)
        if Are_Isom_Mod(Mw,C_mod):
            res.append(w)
    return(res)



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
