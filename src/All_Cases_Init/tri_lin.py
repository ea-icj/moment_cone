import itertools
from sage.all import  identity_matrix

from .root import *
from .weight import *
from .tau import *
from .rep import *



       
def Is_Lin_Triangular(V: Representation,tau:Tau, Invs : list[Root]) -> bool:
    """
    Check if a given list of equalities is Linear Trinagular or not.
    Each equality is actually given by its list of monomials (without coefficients)
    """
    if Invs==[] :
        return(True)
    MG=V.Matrix_Graph(Invs)
    # TODO : modifier avec G^2 là où c'est suffisant ci-dessous
    In= identity_matrix(QQ, V.dim)
    IM=(In-MG).inverse()-In-MG
    listnp=[]
    for ll in list(tau.non_positive_weights(V).values()):
        listnp+=ll    
    Indices_neg_weights=[chi.idx(V) for chi in listnp]  
    cpt=0
    List_inv=[]
    listp=[]
    for ll in list(tau.positive_weights(V).values()):
        listp+=ll
    shifts=[sum(V.G[:i]) for i in range(len(V.G)+1)]
    
    LInv=[]
    for chi in listp:
        id_chi=chi.idx(V)
        if all(IM[id_chi,j]<=1 for j in Indices_neg_weights) : # the equation associated to chi is linear
            cpt+=1
            for j in range(V.dim):
                if IM[id_chi,j] !=0 :
                    # TODO : ce qui suit n'est pas très efficace. A améliorer ?
                    chi2=V.all_weights[j]
                    alpha=chi2.as_vector-chi.as_vector
                    i,j=[k for k in range(len(alpha)) if alpha[k] != 0][:2]
                    k=0
                    while shifts[k]<i: k+=1
                    #alpha=Root(k-1,i-shifts[k-1],j-shifts[k-1])    
                    LInv.append((k-1,i-shifts[k-1],j-shifts[k-1]))
    LInv=list(set(LInv)) 
    LInv=[Root(*l) for l in LInv]
    if len(LInv) != cpt:
        print('tau',tau,Invs,LInv)
    NInv=[alpha for alpha in Invs if alpha not in LInv]
    return(Is_Lin_Triangular(V,tau, NInv))
