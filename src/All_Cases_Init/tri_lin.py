import itertools

from .root import *
from .weight import *
from .tau import *
from .rep import *



       
def Is_Lin_Triangular(V: Representation,tau:Tau, Invs : list[Root]) -> bool:
    """
    Check if a given list of equalities is Linear Trinagular or not.
    Each equality is actually given by its list of monomials (without coefficients)
    """
    if Invs=[] :
        return(True)
    MG=V.Matrix_Graph(Invs)
    # TODO : modifier avec G^2 là où c'est suffisant ci-dessous
    In=M = identity_matrix(QQ, V.dim)
    IM=(In-MG).inverse()-In-MG
    Indeces_neg_weights=[id_chi=chi.idx(V) for chi in tau.negative_weights(V)]
    cpt=0
    List_inv=[]
    Lchi=[chi for chi in tau.positive_weights(V)] # TODO : this is used to get a weight from its index. On doit pouvoir faire mieux
    shifts=[sum(V.G[:i]) for i in range(len(G)+1)]
    #### end TODO
    LInvs=[]
    for chi in tau.positive_weights(V):
        id_chi=chi.idx(V)
        if all(IM[id_chi,j]<=1 for j in Indeces_neg_weights) : # the equation associated to chi is linear
            cpt+=1
            for j in range(V.dim):
                if IM[id_chi,j] !=0 :
                    # TODO : ce qui suit n'est pas très efficace. A améliorer ?
                    chi2=Lchi[j]
                    alpha=chi2.as_vector-chi.as_vector
                    i,j=[k for k in range(len(alpha)) if v[k] != 0][:2]
                    k=0
                    while shifts[k]<i: k+=1
                    #alpha=Root(k-1,i-shifts[k-1],j-shifts[k-1])    
                    LInv.append((k-1,i-shifts[k-1],j-shifts[k-1]))
    LInv=list(set(LInv))
    LInv=[Root(l) for l in LInv]
    if len(LInv) != cpt:
        print('tau',tau,Invs,LInv)
    NInv=[alpha for alpha in Invs if alpha not in LInv]
    return(Is_Lin_Triangular(V,tau, NInvs)
