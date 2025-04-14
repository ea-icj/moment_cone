__all__ = (
    'is_not_contracted',
    'Compute_JA_square_free',
    'Is_Ram_contracted',
)

from random import randint
import itertools
import numpy as np
from numpy.typing import NDArray
import sympy as sp

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
from .utils import prod,fl_dic

    
def is_not_contracted(
        inversions_v: Sequence[Root],
        V: Representation,
        method: Method,
        non_positive_weights: Sequence[Weight],
        positive_weights: Sequence[Weight],
    ) -> bool:
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
    invs_idx=[a.index_in_all_of_U(V.G) for a in inversions_v]
    npw_idx=[V.index_of_weight(chi) for chi in non_positive_weights]
    pw_idx=[V.index_of_weight(chi) for chi in positive_weights]
    #v = point_vect(non_positive_weights, V, ring, bounds=(-10, 10))
    #A = matrix(ring, len(positive_weights), len(inversions_v))
    #for j, root in enumerate(inversions_v):
    #    uv = V.action_op_el(root, v)
    #    for i, chi in enumerate(positive_weights):
    #        A[i, j] = uv[V.index_of_weight(chi)]
    An: NDArray[Any]
    rank_A: int

    if method == "probabilistic" :
        An = V.T_Pi_3D(method, "imaginary")[np.ix_([0,1],npw_idx, pw_idx, invs_idx)].sum(axis=1)
        # Conversion in a sympy matrix of rationnals
        A = sp.Matrix( # type: ignore
            len(pw_idx),
            len(invs_idx),
            lambda i, j:
                    sp.Rational(An[0,i, j]) + sp.Rational(An[1,i, j]) * sp.I
        )       
        rank_A = A.rank() # type: ignore
    else :
        An = V.T_Pi_3D(method, "imaginary")[np.ix_(npw_idx, pw_idx, invs_idx)].sum(axis=0)    
        # Sage matrix
        A = matrix(ring,An)
        rank_A = A.rank() # type: ignore
    
    #A=matrix(ring,An)     
    #rank_A: int = A.change_ring(ring.fraction_field()).rank()
    return rank_A == len(invs_idx)


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
    #v = point_vect(zero_weights, V, ring, bounds=(-10, 10)) # bounds unuseful here
    zw_idx=[V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)]
    gr = ineq.gr_inversions
    Jred: Polynomial = 1
    J: Polynomial = 1
    factors_Jred: list[Polynomial] = []
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau
        gr_idx=[a.index_in_all_of_U(V.G) for a in gr[x]]
        gw_idx=[V.index_of_weight(chi) for chi in tau.positive_weights(V)[x]]
        Mn = V.T_Pi_3D('symbolic')[np.ix_(zw_idx, gw_idx, gr_idx)].sum(axis=0)
        M=matrix(ring,Mn)
        #for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
        #    uv=V.action_op_el(root, v)
        #    for row, chi in enumerate(tau.positive_weights(V)[x]): # List of weights such that tau.scalar(chi)=x 
        #        M[row,col]=uv[V.index_of_weight(chi)]
        
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
    
    gr_rootU = tau.grading_rootsU
    sorted_weightsU = sorted(gr_rootU.keys())
    fl_inv_w=fl_dic(ineq.gr_inversions,sorted_weightsU)

    #print('entry in boundary')
    ### Divisors of the boudary
    for k,w in enumerate(ws):
        for v in w.covering_relations_strong_Bruhat:
            if v.is_min_rep(tau.reduced.mult[k]): 
                vs = list(ws[:k]) + [v] + list(ws[k+1:])
                ineqv = Inequality(tau,w=tuple(vs))
                fl_inv_v=fl_dic(ineqv.gr_inversions,sorted_weightsU)
                if any([fl_inv_v[key] > fl_inv_w[key] for key in fl_inv_v.keys()]):
                    continue
                #for i,j in v.inversions : 
                #    if Root(k,i,j) not in ineq.inversions : ## TODO : c'est all(alha + chi not in weights for chi orthogonal)
                #    h=tau.dot_root(Root(k,i,j))
                #    if [fl_inv_v[h-1] == fl_inv_w[h-1]:
       	        #        continue
                if is_not_contracted(tuple(ineqv.inversions),V,method_S,Neg0_Weights_sorted,Pos_Weights_sorted) :
                    return(False)

    ### Divisor R_0
    J,J_square_free, factors_J_sqf= Compute_JA_square_free(ineq, V) # The Jacobian and it's reduced form
    
    # Generic point v of V(tau<=0) and matrix of Tpi at (e,v)
    #v = point_vect(Neg0_Weights_sorted,V,V.QV)
    #A=matrix(V.QV,len(Pos_Weights_sorted),dU)
    #B0=matrix(V.QV,len(tau.orthogonal_weights(V)),dU)
    #gr = ineq.gr_inversions

    #col=0
    #for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau
    #    if x > 0 :
    #        gr_idx=[a.index_in_all_of_U(G) for a in gr[x]]
    #        gw_idx=[V.index_of_weight(chi) for chi in gw[x]]   
    #        for root in gr[x]: # List of roots such that tau.scalar(root)=x
    #            uv=V.action_op_el(root, v)
    #            for row,chi in enumerate(Pos_Weights_sorted):
    #                A[row,col]=uv[V.index_of_weight(chi)]
    #            for row,chi in enumerate(tau.orthogonal_weights(V)):
    #                B0[row,col]=uv[V.index_of_weight(chi)]
    #            col+=1
       
    inv_idx=[a.index_in_all_of_U(V.G) for a in ineq.inversions]
    pw_idx=[V.index_of_weight(chi) for chi in Pos_Weights_sorted]
    npw_idx=[V.index_of_weight(chi) for chi in Neg0_Weights_sorted]
    zw_idx=[V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)]
    Azn = V.T_Pi_3D(method_R0, 'line')[np.ix_(npw_idx, pw_idx, inv_idx)].sum(axis=0)
    Az=matrix(ring_R0,Azn)
    B0zn = V.T_Pi_3D(method_R0, 'line')[np.ix_(npw_idx, zw_idx, inv_idx)].sum(axis=0) 
    B0z=matrix(ring_R0,B0zn)
    # The line: gradiant of J
    L0=matrix(V.QV,1,len(tau.orthogonal_weights(V)))
    for col,idx in enumerate(zw_idx):
        chi = V.all_weights[idx]
        L0[0,col]=J_square_free.derivative(V.QV.variable(chi))
    #print('L0 done')    
    #for col,chi in enumerate(tau.orthogonal_weights(V)) :
    #    L0[0,col]=J_square_free.derivative(V.QV.variable(chi))

    # Dictionnary for substitution    
    #subs_dict: dict[Variable, Variable] = {}    
    #for chi in Neg0_Weights_sorted:
    #    if method_R0 == "probabilistic":
    #        subs_dict[V.QV.variable(chi)]= randint(-500,500)*V.QZ('z')+randint(-500,500)# TODO :Tester l'effet du changement de 500. Math : Doit-on mettre du I ? 
    #    else:
    #        va, vb = V.QV2.variable(chi) 
    #        subs_dict[V.QV.variable(chi)]= va*ring_R0('z') + vb # type: ignore

    # Substitutions
    #Az=A.subs(subs_dict)
    
    #B0z=B0.subs(subs_dict)
    L0z=L0.subs(V.T_Pi_3D(method_R0, 'dict'))
    #Jz=J.subs(subs_dict)
    Jz=J.subs(V.T_Pi_3D(method_R0, 'dict'))
    assert method_R0 == 'symbolic' or J.degree() == Jz.degree(), "The random line is not enough generic to intersect each irreducible component of R0. Please Restart."
    #print("L0z Jz done")
    factors_J_sqf_z = sum([list(dict(Poly.subs(V.T_Pi_3D(method_R0, 'dict')).factor()).keys()) for Poly in factors_J_sqf],[])
    #print("factor J done")
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
    #print("delta done")
    # Computation of Bezout inverse
    LIB=Bezout_Inverse(Ldelta,ring_R0)
    #print("Bezout done")
    # Kernel of Az modulo delta
    noyau=Kernel_modulo_P(ring_R0,Az,Ldelta,LIB)
    #print("kernel done")
    #print("inversions",ineq.inversions)
    #print("orth weights",tau.orthogonal_weights(V))
    #print("positive weights",tau.positive_weights(V))
    #if Inequality.from_tau(ineq.wtau.end0_representative.sort_mod_sym_dim)  == Inequality.from_tau(Tau(((-1, 0, 1, 0), (0, 1, -1, 0), (1, 0, -1, 0), (-1,))).end0_representative.sort_mod_sym_dim):
    #    print('kernel',noyau)
    # Check divisibility
    check=L0z*B0z*noyau
    check=ring_R0(check[0])
    quo, rem = check.quo_rem(delta)
    return cast(bool, rem == 0)
   

