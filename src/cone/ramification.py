from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .cone_dimension import point_vect, vector, matrix
from .spaces import action_op_el
from .inequality import *
from .permutation import *

# FIXME: we get d from tau but in the current code, it will leads to recreate the rings for each tau.
def is_not_contracted(inverse_v: Sequence[Root], tau: Tau, method: Method) -> bool:
    """
    ???

    The vector v is chosen depending on the method: probabilistic or symbolic.

    FIXME: example

    Example:
    >>> from cone import *
    >>> d = Dimension((2, 3, 4))
    >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
    >>> tau
    1 | 6 2 | 1 4 1 | 4 5 3 1
    >>> roots = list(Root.all_of_U(d))[:10]
    >>> is_not_contracted(roots, tau, "probabilistic")
    False
    >>> is_not_contracted(roots, tau, "symbolic")
    False
    """
    d = tau.d

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = d.QI
    elif method == "symbolic":
        ring = d.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    
    # FIXME: do we stay we list conversion at each call?
    # Maybe grading root and weight should be implemented using a more convenient class?
    from itertools import chain
    non_positive_weights = list(chain.from_iterable(tau.non_positive_weights.values()))
    positive_weights = list(chain.from_iterable(tau.positive_weights.values()))

    # FIXME: why not positive_weights instead of non_negative_weights?
    v = point_vect(non_positive_weights, d, ring, bounds=(-1000, 1000))
    list_inverse_v=list(inverse_v)
    T = matrix(ring, len(positive_weights), len(list_inverse_v))
    for j, root in enumerate(list_inverse_v):
        uv = action_op_el(root, v, d)
        for i, chi in enumerate(positive_weights):
            T[i, j] = uv[chi.index_in(d)]

    rank_T = T.change_ring(ring.fraction_field()).rank()
    return rank_T == len(list_inverse_v)

def Compute_JA(ineq : Inequality) : # Return a dictionnary polynom :  int
    tau=ineq.tau
    d = tau.d
    ring = d.QV
    # a generic vector in VV^tau
    zero_weights = tau.orthogonal_weights
    v = point_vect(zero_weights, d, ring, bounds=(-100, 100)) # bounds unuseful here
    # inversions of w
    Inv_w=ineq.inversions
    #gr = grading_dictionary(ineq.inversions, tau.dot_root)
    gr = tau.grading_roots_in(ineq.inversions)
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau        
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
            uv=action_op_el(root, v, d)
            for row, chi in enumerate(gw[x]): # List of weights such that tau.scalar(chi)=x 
                M[row,col]=uv[chi.index_ind(d)]
        Jb=M.det().factor()   
        for F in Jb.keys(): # We could make a function add_dictionaries
            if F in J.keys():
                J[F]+=Jb[F]
            else:
                J[F]=Jb[F]
    return(J)

def Smith_n_1(A):
    "Compute the gcd of the minors of A of size n-1"
    combinaisons = list(combinations(range(A.nrows()),A.nrows()-1))
    pgcd = 0
    
    # Run over the n-1 x n-1 minors
    for lignes in combinaisons:
        for colonnes in combinaisons:
            # Extract submatrices
            sous_matrice = A.matrix_from_rows_and_columns(lignes, colonnes)
            
            # Det            
            det = sous_matrice.det()
            
            # Upgrade GCD
            pgcd = gcd(pgcd, det)
            
            #If the gcd is 1 stop
            if pgcd == 1:
                return pgcd    
    return pgcd

def Is_Ram_contracted(ineq : Inequality, method_S: Method, method_R0: Method) -> bool :
    d=ineq.tau.d
    if method_S == "probabilistic":
        ring = d.QI
    elif method_S == "symbolic":
        ring = d.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    ws=ineq.w
    Inv_w=list(ineq.inversions)
    dU=len(Inv_w)
    tau=ineq.tau
    
    # Creation of sorted lists of weights
    Neg0_Weights_dic=tau.non_positive_weights
    Neg0_Weights_sorted=[]
    for x in sorted(Neg0_Weights_dic.keys(),reverse=True):
        Neg0_Weights_sorted+=Neg0_Weights_dic[x]
    Pos_Weights_dic=tau.positive_weights
    Pos_Weights_sorted=[]
    for x in sorted(Pos_Weights_dic.keys(),reverse=True):
        Pos_Weights_sorted+=Pos_Weights_dic[x]
    if dU==0 : # Trivial case. Can we replace by <=1
        return(True)

    ### Divisors of the boudary
    for k,w in enumerate(ws):
        for v in w.covering_relations_strong_Bruhat:
            if v.is_min_rep(tau.reduced.mult[k]): 
                vs=Permutation(ws[:k]+(v,)+ws[k+1:])
                ineqv=Inequality(tau,vs)
                if is_not_contracted(ineqv.inversions,tau,method_S) :
                    return(False)
    return(True) #### A SUPPRIMER QUAND ON AURA R0#####
