from random import randint
import itertools
from sage.all import PolynomialRing

from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .rep import *
from .vector_chooser import point_vect, vector, matrix
from .inequality import *
from .permutation import *
from .kx_mod import *

__all__ = (
    'is_not_contracted',
    'Compute_JA',
    'Is_Ram_contracted',
)

    
# FIXME: we get d from tau but in the current code, it will leads to recreate the rings for each tau.
def is_not_contracted(inversions_v: Iterable[Root], tau: Tau, V: Representation,method: Method) -> bool:
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

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = V.QI
    elif method == "symbolic":
        ring = V.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    
    # FIXME: do we stay we list conversion at each call?
    # Maybe grading root and weight should be implemented using a more convenient class?
    from itertools import chain
    non_positive_weights = list(chain.from_iterable(tau.non_positive_weights(V).values())) # todo : Cela est aussi fait avant l'appel à la fonction.
    positive_weights = list(chain.from_iterable(tau.positive_weights(V).values())) # todo : Cela est aussi fait avant l'appel à la fonction.

    v = point_vect(non_positive_weights, V, ring, bounds=(-1000, 1000))
    list_inversions_v = list(inversions_v)
    A = matrix(ring, len(positive_weights), len(list_inversions_v))
    for j, root in enumerate(list_inversions_v):
        uv = V.action_op_el(root, v)
        for i, chi in enumerate(positive_weights):
            A[i, j] = uv[chi.idx(V)]

    rank_A = A.change_ring(ring.fraction_field()).rank()
    return rank_A == len(list_inversions_v)

def Normalization_Factorized_Polynomial(Jb) :
    d={}
    for P in Jb.keys():
        a=P.monomial_coefficient(P.monomials()[0])
        new_key=P/a
        d[new_key]=Jb[P]
    return(d)

def Compute_JA(ineq : Inequality,V: Representation) : # Return a dictionnary polynom :  int
    tau=ineq.tau
    ring = V.QV
    # a generic vector in VV^tau
    zero_weights = tau.orthogonal_weights(V)
    v = point_vect(zero_weights, V, ring, bounds=(-100, 100)) # bounds unuseful here
    # inversions of w
    Inv_w=ineq.inversions
    #gr = grading_dictionary(ineq.inversions, tau.dot_root)
    gr = tau.grading_roots_in(ineq.inversions)
    J: dict[Any, Any] = {} # FIXME type
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau        
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
            uv=V.action_op_el(root, v)
            for row, chi in enumerate(tau.positive_weights(V)[x]): # List of weights such that tau.scalar(chi)=x 
                M[row,col]=uv[chi.idx(V)]
    
        Jb=dict(M.det().factor())
        Jbn=Normalization_Factorized_Polynomial(Jb)
        for F in Jbn.keys(): # We could make a function add_dictionaries
            if F in J.keys():
                J[F]+=Jbn[F]
            else:
                J[F]=Jbn[F]
    return(J)

def Smith_n_1(A):
    "Compute the gcd of the minors of A of size n-1"
    combinaisons = list(itertools.combinations(range(A.nrows()),A.nrows()-1))
    ring=A.base_ring()
    pgcd = ring(0)
    
    # Run over the n-1 x n-1 minors
    for lignes in combinaisons:
        for colonnes in combinaisons:
            # Extract submatrices
            sous_matrice = A.matrix_from_rows_and_columns(lignes, colonnes)
            
            # Det            
            det_Smith = sous_matrice.det()
            
            # Upgrade GCD
            pgcd = pgcd.gcd(det_Smith)
            
            #If the gcd is 1 stop
            if pgcd == 1:
                return pgcd    
    return pgcd

def Is_Ram_contracted(ineq : Inequality, V: Representation, method_S: Method, method_R0: Method) -> bool :
    ws=ineq.w
    Inv_w=list(root for root in ineq.inversions)
    dU=len(Inv_w)
    if dU==0 : # Trivial case. TODO : Math : Can we replace by <=1 ?
        return(True)
    
    if method_R0 == "probabilistic" :
        ring_R0= V.QZ
    elif method_R0 == "symbolic":
        K=V.QV2.fraction_field()
        ring_R0 = PolynomialRing(K,"z")
    else:
        raise ValueError(f"Invalid value {method_R0} of the computation method")

    tau=ineq.tau
    
    # Creation of sorted lists of weights todo: modifier en utilisant mieux la classe tau
    # TODO : ici on ordonne les clés et crée une liste. Peut-être on peut faire plus simple et dans utils
    Neg0_Weights_dic=tau.non_positive_weights(V) 
    Neg0_Weights_sorted=[]
    for x in sorted(Neg0_Weights_dic.keys(),reverse=True):
        Neg0_Weights_sorted+=Neg0_Weights_dic[x]
    Pos_Weights_dic=tau.positive_weights(V)
    Pos_Weights_sorted=[]
    for x in sorted(Pos_Weights_dic.keys(),reverse=True):
        Pos_Weights_sorted+=Pos_Weights_dic[x]
    

    ### Divisors of the boudary
    for k,w in enumerate(ws):
        for v in w.covering_relations_strong_Bruhat:
            if v.is_min_rep(tau.reduced.mult[k]): 
                vs = list(ws[:k]) + [v] + list(ws[k+1:])
                ineqv = Inequality(tau,vs)
                if is_not_contracted(ineqv.inversions,tau,V,method_S) :
                    return(False)

    ### Divisor R_0
    Jf=Compute_JA(ineq,V) # The Jacobian factorized as a dictionnary
    
    J_square_free=1
    for pol in Jf.keys():
        J_square_free*=pol # todo : prod(list(Jf.keys())) ne semble pas fonctionner
    # FIXME: type ignore
    if len(Jf.keys())!=len(dict(J_square_free.factor()).keys()): # type: ignore
        print('Error in factor with:',Jf,J_square_free)
 
    # Generic point v of V(tau<=0) and matrix of Tpi at (e,v)
    v = point_vect(Neg0_Weights_sorted,V,V.QV)
    A=matrix(V.QV,len(Pos_Weights_sorted),dU)
    B0=matrix(V.QV,len(tau.orthogonal_weights(V)),dU)
    gr = tau.grading_roots_in(ineq.inversions)
    col=0
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau
        if x > 0 :   
            for root in gr[x]: # List of roots such that tau.scalar(root)=x
                uv=V.action_op_el(root, v)
                for row,chi in enumerate(Pos_Weights_sorted):
                    A[row,col]=uv[chi.idx(V)]
                for row,chi in enumerate(tau.orthogonal_weights(V)):
                    B0[row,col]=uv[chi.idx(V)]
                col+=1
       
    # The line: gradiant of J
    L0=matrix(V.QV,1,len(tau.orthogonal_weights(V)))
    for col,chi in enumerate(tau.orthogonal_weights(V)) :
        L0[0,col]=J_square_free.derivative(V.QV.variable(chi))

    # Dictionnary for substitution    
    subs_dict = {}    
    for chi in Neg0_Weights_sorted:
        if method_R0 == "probabilistic":
            subs_dict[V.QV.variable(chi)]= randint(-500,500)*V.QZ('z')+randint(-500,500)# TODO :Tester l'effet du changement de 500. Math : Doit-on mettre du I ? 
        else:
            va, vb = V.QV2.variable(chi) 
            subs_dict[V.QV.variable(chi)]= va*ring_R0('z')+vb 
    # Substitutions
    Az=A.subs(subs_dict)
   
    B0z=B0.subs(subs_dict)
    L0z=L0.subs(subs_dict)
    # FIXME: type ignore
    Jz=J_square_free.subs(subs_dict) # type: ignore


    # Computation of reduced delta as factorized polynomial
    Smith_n_un=Smith_n_1(Az)
    # FIXME: type ignore
    s=Jz.degree()-Smith_n_un.degree() # type: ignore
    exponent=max(s-1,1)
    pgcd=Jz.gcd(Smith_n_un**exponent)
    
    delta1 = Jz // pgcd 
    delta = delta1 // delta1.gcd(delta1.derivative()) 
    Ldelta=list(dict(delta.factor()).keys())
    

    # Computation of Bezout inverse
    LIB=Bezout_Inverse(Ldelta,ring_R0)
    
    # Kernel of Az modulo delta
    noyau=Kernel_modulo_P(ring_R0,Az,Ldelta,LIB)

    # Check divisibility
    check=L0z*B0z*noyau
    check=ring_R0(check[0])
    quo, rem = check.quo_rem(delta)
    return(rem==0)
   

