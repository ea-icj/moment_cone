


__all__ = (
    'is_not_contracted',
    'Compute_JA_square_free',
    'Is_Ram_contracted',
)

from random import randint
import itertools

from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .representation import *
from .vector_chooser import point_vect, vector, matrix, Matrix
from .inequality import *
from .permutation import *
from .kx_mod import *
from .rings import PolynomialRing, Polynomial, Variable
from .utils import prod

    
# FIXME: we get d from tau but in the current code, it will leads to recreate the rings for each tau.
def is_not_contracted(inversions_v: Iterable[Root], tau: Tau, V: Representation, method: Method) -> bool:
    """
    ???

    The vector v is chosen depending on the method: probabilistic or symbolic.

    FIXME: example

    Example:
    >>> from cone import *
    >>> G = LinearGroup((2, 3, 4, 1))
    >>> V = KroneckerRepresentation(G)
    >>> tau = Tau.from_flatten([6, 2, 1, 4, 1, 4, 5, 3, 1, 1], G)
    >>> tau
    6 2 | 1 4 1 | 4 5 3 1 | 1
    >>> roots = list(Root.all_of_U(G))[:10]
    >>> is_not_contracted(roots, tau, V, "probabilistic")
    False
    >>> is_not_contracted(roots, tau, V, "symbolic")
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
            A[i, j] = uv[V.index_of_weight(chi)]
    #print('inv',list_inversions_v)
    #print('pos weights', positive_weights)
    #print('A',A)        

    rank_A: int = A.change_ring(ring.fraction_field()).rank()
    return rank_A == len(list_inversions_v)

def Normalization_Factorized_Polynomial(Jb: dict[Polynomial, int]) -> dict[Polynomial, int]:
    d: dict[Polynomial, int] = {}
    for P in Jb.keys():
        a: int = P.monomial_coefficient(P.monomials()[0])
        new_key=P/a
        d[new_key]=Jb[P]
    return d

def Compute_JA_square_free(ineq: Inequality, V: Representation) -> (Polynomial,Polynomial,Iterable[Polynomial]):
    tau = ineq.tau
    ring = V.QV
    # a generic vector in VV^tau
    zero_weights = tau.orthogonal_weights(V)
    v = point_vect(zero_weights, V, ring, bounds=(-100, 100)) # bounds unuseful here
    #gr = grading_dictionary(ineq.inversions, tau.dot_root)
    gr = tau.grading_roots_in(ineq.inversions)
    Jred: Polynomial = 1
    J: Polynomial = 1
    factors_Jred=[]
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau        
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
            uv=V.action_op_el(root, v)
            for row, chi in enumerate(tau.positive_weights(V)[x]): # List of weights such that tau.scalar(chi)=x 
                M[row,col]=uv[V.index_of_weight(chi)]
        #print('M',M)
        
        Jb: Polynomial = M.det()
        partial_derivatives: list[Polynomial] = [
            Jb.derivative(V.QV.variable(chi))
            for chi in zero_weights
        ]
        from sage.all import gcd as sage_gcd # type: ignore
        Jbred_denom: Polynomial = sage_gcd([Jb] + partial_derivatives)
        Jbred: Polynomial
        remainder: Polynomial
        Jbred, remainder = Jb.quo_rem(Jbred_denom)
        assert remainder == 0
        
        Jred_new, remainder = (Jred * Jbred).quo_rem(sage_gcd(Jred, Jbred)) #computing the square free version of Jred*Jbred
        factor_Jred,remainder2=Jred_new.quo_rem(Jred)
        if factor_Jred.degree()!=0:
            factors_Jred.append(factor_Jred) #new factor in Jred at this level
        J *= Jb
        Jred=Jred_new
        assert (remainder == 0) and (remainder2==0)

    return J,Jred,factors_Jred

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
                #print('Schub Div',vs)
                ineqv = Inequality(tau,vs)
                if is_not_contracted(ineqv.inversions,tau,V,method_S) :
                    return(False)

    ### Divisor R_0
    J,J_square_free, factors_J_sqf= Compute_JA_square_free(ineq, V) # The Jacobian and it's reduced form
    
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
                    A[row,col]=uv[V.index_of_weight(chi)]
                for row,chi in enumerate(tau.orthogonal_weights(V)):
                    B0[row,col]=uv[V.index_of_weight(chi)]
                col+=1
       
    # The line: gradiant of J
    L0=matrix(V.QV,1,len(tau.orthogonal_weights(V)))
    for col,chi in enumerate(tau.orthogonal_weights(V)) :
        L0[0,col]=J_square_free.derivative(V.QV.variable(chi))

    # Dictionnary for substitution    
    subs_dict: dict[Variable, Variable] = {}    
    for chi in Neg0_Weights_sorted:
        if method_R0 == "probabilistic":
            subs_dict[V.QV.variable(chi)]= randint(-500,500)*V.QZ('z')+randint(-500,500)# TODO :Tester l'effet du changement de 500. Math : Doit-on mettre du I ? 
        else:
            va, vb = V.QV2.variable(chi) 
            subs_dict[V.QV.variable(chi)]= va*ring_R0('z')+vb # type: ignore

    # Substitutions
    Az=A.subs(subs_dict)
   
    B0z=B0.subs(subs_dict)
    L0z=L0.subs(subs_dict)
    Jz=J.subs(subs_dict)
    factors_J_sqf_z = sum([list(dict(Poly.subs(subs_dict).factor()).keys()) for Poly in factors_J_sqf],[])
    #print(factors_J_sqf_z)
    #print(list(dict(J_sqf_z.factor()).keys()))
    #J_sqf_z=J_square_free.subs(subs_dict)
    #factors_J_sqf_z=list(dict(J_sqf_z.factor()).keys())
    Ldelta=[]
    for delta1 in factors_J_sqf_z:
        quo, rem = Jz.quo_rem(delta1**2)
        if rem != 0 :
            Ldelta.append(delta1)
        else :
            delta1_quotient = ring_R0.quotient(delta1)
            Ared=Az.apply_map(lambda entry: delta1_quotient(entry)) # A modulo delta1
            if Ared.rank() == Ared.ncols()-1 :
                Ldelta.append(delta1)
    if Ldelta==[]:
        delta=1
    else: 
        delta = prod(Ldelta)            
    

    # Computation of Bezout inverse
    LIB=Bezout_Inverse(Ldelta,ring_R0)
    
    # Kernel of Az modulo delta
    noyau=Kernel_modulo_P(ring_R0,Az,Ldelta,LIB)

    # Check divisibility
    check=L0z*B0z*noyau
    check=ring_R0(check[0])
    quo, rem = check.quo_rem(delta)
    return cast(bool, rem == 0)
   

