


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

    
def is_not_contracted(inversions_v: tuple[Root], tau: Tau, V: Representation, method: Method,non_positive_weights: list[Weight],positive_weights: list[Weight]) -> bool:

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

    v = point_vect(non_positive_weights, V, ring, bounds=(-10, 10))
    #list_inversions_v = list(inversions_v)
    A = matrix(ring, len(positive_weights), len(inversions_v))
    for j, root in enumerate(inversions_v):
        uv = V.action_op_el(root, v)
        for i, chi in enumerate(positive_weights):
            A[i, j] = uv[V.index_of_weight(chi)]   

    #print(A)        

    rank_A: int = A.change_ring(ring.fraction_field()).rank()
    #print(rank_A)
    return rank_A == len(inversions_v)

def Normalization_Factorized_Polynomial(Jb: dict[Polynomial, int]) -> dict[Polynomial, int]:
    d: dict[Polynomial, int] = {}
    for P in Jb.keys():
        a: int = P.monomial_coefficient(P.monomials()[0])
        new_key=P/a
        d[new_key]=Jb[P]
    return d

def Compute_JA_square_free(ineq: Inequality, V: Representation) -> tuple[Polynomial, Polynomial, list[Polynomial]]:
    tau = ineq.tau
    ring = V.QV
    # a generic vector in VV^tau
    zero_weights = tau.orthogonal_weights(V)
    v = point_vect(zero_weights, V, ring, bounds=(-10, 10)) # bounds unuseful here
    gr = ineq.gr_inversions
    Jred: Polynomial = 1
    J: Polynomial = 1
    factors_Jred: list[Polynomial] = []
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau        
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
            uv=V.action_op_el(root, v)
            for row, chi in enumerate(tau.positive_weights(V)[x]): # List of weights such that tau.scalar(chi)=x 
                M[row,col]=uv[V.index_of_weight(chi)]
        
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

    return J, Jred, factors_Jred

def Is_Ram_contracted(ineq : Inequality, V: Representation, method_S: Method, method_R0: Method) -> bool :
    ineq_check=Inequality.from_tau(Tau(((0, -1, 0, 1),(0, 1, -1, 0),(1, 0, -1, 0),(-1,))))
    #if Inequality.from_tau(ineq.wtau.end0_representative.sort_mod_sym_dim)  == Inequality.from_tau(Tau(((-1, 0, 1, 0), (0, 1, -1, 0), (1, 0, -1, 0), (-1,))).end0_representative.sort_mod_sym_dim):
    #    print('coucou Ram',ineq, ineq.inversions, ineq.gr_inversions)
    #print('Tau',ineq.tau)    
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
    
    # Creation of sorted lists of weights 
    #Neg0_Weights_dic=tau.non_positive_weights(V) 
    Neg0_Weights_sorted=list(itertools.chain.from_iterable(tau.non_positive_weights(V)[k] for k in sorted(tau.non_positive_weights(V))))
    Pos_Weights_sorted=list(itertools.chain.from_iterable(tau.positive_weights(V)[k] for k in sorted(tau.positive_weights(V))))
    
    #print('entry in boundary')
    ### Divisors of the boudary
    for k,w in enumerate(ws):
        for v in w.covering_relations_strong_Bruhat:
            if v.is_min_rep(tau.reduced.mult[k]): 
                vs = list(ws[:k]) + [v] + list(ws[k+1:])
                #print('Schub Div',vs)
                #print('v',vs[0])
                #print(vs[0].inversions)
                #print(type(vs[0]))
                ineqv = Inequality(tau,w=tuple(vs))
                #print(ineq)
                #print(ineqv)
                #print(ineqv.inversions)
                #print(len(Neg0_Weights_sorted))
                #print(len(Pos_Weights_sorted))
                if is_not_contracted(ineqv.inversions,tau,V,method_S,Neg0_Weights_sorted,Pos_Weights_sorted) :
                    return(False)
                
    #if Inequality.from_tau(ineq.wtau.end0_representative.sort_mod_sym_dim)  == Inequality.from_tau(Tau(((-1, 0, 1, 0), (0, 1, -1, 0), (1, 0, -1, 0), (-1,))).end0_representative.sort_mod_sym_dim):
    #    print('Schub passed')            
    #print('boundary passed')
    ### Divisor R_0
    J,J_square_free, factors_J_sqf= Compute_JA_square_free(ineq, V) # The Jacobian and it's reduced form
    
    # Generic point v of V(tau<=0) and matrix of Tpi at (e,v)
    v = point_vect(Neg0_Weights_sorted,V,V.QV)
    A=matrix(V.QV,len(Pos_Weights_sorted),dU)
    B0=matrix(V.QV,len(tau.orthogonal_weights(V)),dU)
    gr = ineq.gr_inversions
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
    
    #if Inequality.from_tau(ineq.wtau.end0_representative.sort_mod_sym_dim)  == Inequality.from_tau(Tau(((-1, 0, 1, 0), (0, 1, -1, 0), (1, 0, -1, 0), (-1,))).end0_representative.sort_mod_sym_dim):
    #    print('delta',delta)

    # Computation of Bezout inverse
    LIB=Bezout_Inverse(Ldelta,ring_R0)
    
    # Kernel of Az modulo delta
    noyau=Kernel_modulo_P(ring_R0,Az,Ldelta,LIB)

    #if Inequality.from_tau(ineq.wtau.end0_representative.sort_mod_sym_dim)  == Inequality.from_tau(Tau(((-1, 0, 1, 0), (0, 1, -1, 0), (1, 0, -1, 0), (-1,))).end0_representative.sort_mod_sym_dim):
    #    print('kernel',noyau)
    # Check divisibility
    check=L0z*B0z*noyau
    #L0B0z=(L0*B0).subs(subs_dict)
    #check=L0B0z*noyau
    check=ring_R0(check[0])
    quo, rem = check.quo_rem(delta)
    #if Inequality.from_tau(ineq.wtau.end0_representative.sort_mod_sym_dim)  == Inequality.from_tau(Tau(((-1, 0, 1, 0), (0, 1, -1, 0), (1, 0, -1, 0), (-1,))).end0_representative.sort_mod_sym_dim):
    #    print('Ram0',rem)
    #    print(ineq)
    return cast(bool, rem == 0)
   

