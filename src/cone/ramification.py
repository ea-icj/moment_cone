__all__ = (
    'is_not_contracted',
    'Compute_JA_square_free',
    'Is_Ram_contracted',
)

from random import randint
import itertools
import numpy as np
from numpy.typing import NDArray
from sage.all import NumberField, fast_callable
from sage.misc.fpickle import pickle_function

from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .representation import *
from .inequality import *
from .permutation import *
#from .kx_mod import *
from .rings import QQ, I, matrix, vector, Matrix, Polynomial, PolynomialRing
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
    
    invs_idx=[a.index_in_all_of_U(V.G) for a in inversions_v]
    npw_idx=[V.index_of_weight(chi) for chi in non_positive_weights]
    pw_idx=[V.index_of_weight(chi) for chi in positive_weights]
    rank_A: int

    if method == "probabilistic" :
        for p in range(V.random_deep):
            An = V.T_Pi_3D(method, "imaginary")[np.ix_([2*p, 2*p+1], npw_idx, pw_idx, invs_idx)].sum(axis=1)
            A = matrix(
                len(pw_idx),
                len(invs_idx),
                lambda i, j: QQ(An[0, i, j]) + I * QQ(An[1, i, j])
                )
            rank_A = A.rank()
            if rank_A == len(invs_idx):
                    break

    else :
        An = V.T_Pi_3D(method, "imaginary")[np.ix_(npw_idx, pw_idx, invs_idx)].sum(axis=0)    
        # Sage matrix
        A = matrix(ring, An)
        rank_A = A.rank()
    
    return rank_A == len(invs_idx)


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
    gr_inv = ineq.gr_inversions 
    Neg0_Weights_sorted=list(itertools.chain.from_iterable(tau.non_positive_weights(V)[k] for k in sorted(tau.non_positive_weights(V)))) 
    gr_rootU = tau.grading_rootsU
    sorted_weightsU = sorted(gr_rootU.keys())
    Pos_Weights_sorted=list(itertools.chain.from_iterable(tau.positive_weights(V)[x] for x in sorted(gr_inv.keys(),reverse=True)))
    fl_inv_w=fl_dic(gr_inv,sorted_weightsU)
    
    ### Divisors of the boudary
    for k,w in enumerate(ws):
        for v in w.covering_relations_strong_Bruhat:
            if v.is_min_rep(tau.reduced.mult[k]): 
                vs = list(ws[:k]) + [v] + list(ws[k+1:])
                ineqv = Inequality(tau,w=tuple(vs))
                fl_inv_v=fl_dic(ineqv.gr_inversions,sorted_weightsU)
                if any([fl_inv_v[key] > fl_inv_w[key] for key in fl_inv_v.keys()]):
                    continue
                if is_not_contracted(tuple(ineqv.inversions),V,method_S,Neg0_Weights_sorted,Pos_Weights_sorted) :
                    return(False)

    ### Divisor R_0
    # Indices of roots and weights    
    InversionSorted=list(itertools.chain.from_iterable(
            gr_inv[x] for x in sorted(gr_inv.keys(),reverse=True)))
    inv_idx=[a.index_in_all_of_U(V.G) for a in InversionSorted]
    pw_idx=[V.index_of_weight(chi) for chi in Pos_Weights_sorted]
    npw_idx=[V.index_of_weight(chi) for chi in Neg0_Weights_sorted]
    zw_idx=[V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)]

    sizeblocks=[0]
    for x in sorted(gr_inv.keys(),reverse=True): 
        sizeblocks.append(sizeblocks[-1]+len(gr_inv[x]))
    

    for p in range(V.random_deep):
        # Construct the matrices Az and B0z 
        if  method_R0 == 'probabilistic':   
            Azn = V.T_Pi_3D(method_R0,'line')[np.ix_([2*p, 2*p+1],npw_idx, pw_idx, inv_idx)].sum(axis=1)
            Az = matrix(ring_R0, Azn[0, :, :] * ring_R0('z') + Azn[1, :, :])
            B0zn = V.T_Pi_3D(method_R0, 'line')[np.ix_([2*p, 2*p+1], npw_idx, zw_idx, inv_idx)].sum(axis=1)
            B0z = matrix(ring_R0, B0zn[0, :, :] * ring_R0('z') + B0zn[1, :, :])
        else :
            Azn = V.T_Pi_3D(method_R0, 'line')[np.ix_(npw_idx, pw_idx, inv_idx)].sum(axis=0)
            Az=matrix(ring_R0,Azn)
            B0zn = V.T_Pi_3D(method_R0, 'line')[np.ix_(npw_idx, zw_idx, inv_idx)].sum(axis=0) 
            B0z=matrix(ring_R0,B0zn)

        List_C_mult_dec=sorted(gr_inv.keys(),reverse=True)
        # Compute Jz and its factorizations
        Jz=1
        Factors_Jz=[]
        Blocks_Az=[]
        List_deltas=[]
        List_Comatices=[]
        shift=0
        for x in List_C_mult_dec : #sorted(gr_inv.keys(),reverse=True):
            n_block=len(gr_inv[x])
            Azi=Az.submatrix(shift,shift,n_block,n_block)
            Ji=Azi.det()
            Azi_com=Ji*Azi.inverse()
            #print('test com :',Azi.adjugate()==Azi_com)
            #TODO : choisir entre ces deux méthodes
            List_Comatices.append(Azi_com.transpose())
            Jz*=Ji
            Blocks_Az.append(Azi)
            Factors_Jz.append(Ji)
            List_deltas.append(dict(Ji.factor()))
            shift+=n_block
            
        #print(ineq)   
        # Run over the diagonal blocks of A
        List_checked_deltas=[]
        List_C_mult_dec=sorted(gr_inv.keys(), reverse=True)
        for i,x in enumerate(List_C_mult_dec) : #sorted(gr_inv.keys(),reverse=True): 
            # Current block and its determinant
            Azi=Blocks_Az[i]
            Ji=Factors_Jz[i]
            n_current=len(gr_inv[x])
            #deltas=dict(Ji.factor())

            for delta, m in List_deltas[i].items():
                if delta not in List_checked_deltas:     
                    List_checked_deltas.append(delta)
                    # Occurences of delta
                    occ_delta=[[i,m,x]]
                    for j,y in enumerate(List_C_mult_dec[i+1:]): 
                        if delta in List_deltas[j].items() :
                            occ_delta.append([j,List_deltas[j][delta],y])
                    #jmax=max([occ[0] for occ in occ_delta])
                    jmin,mult_min,y = min(occ_delta, key=lambda occ: [occ[1],len(gr_inv[occ[2]])])
                    # We minimize multiplicity of delta. If several ones, we minimize the size of the block                      
                    
                    K=NumberField(delta, 'a')
                    a = K.gen()
                    # submatric of A modulo delta used to compute the kernel
                    Ared=Az.matrix_from_columns(range(sizeblocks[i], Az.ncols())).apply_map(lambda entry: entry(a))
                    B0z_red=B0z.matrix_from_columns(range(sizeblocks[i], Az.ncols())).apply_map(lambda entry: entry(a))
                    # Block used to compoute L0
                    Aredj=Blocks_Az[jmin].apply_map(lambda entry: entry(a)) # A modulo delta
                    com_Aredj=List_Comatices[jmin].apply_map(lambda entry: entry(a)) # A modulo delta
                    
                    if Jz % delta**2 != 0  or Ared.rank() == Ared.ncols()-1 :
                        noyau=Ared.right_kernel().basis()[0] 

                        #noyau_converted = [ring_R0(n.lift()) for n in noyau]
                        # Computation of L0   
                        if mult_min==1:
                            nA = Ared.nrows()
                            L0=matrix(K,1,len(zw_idx))
                            #y=gr_inv.keys()[jmin]
                            
                            for col,idx in enumerate(zw_idx):
                                Mn = V.T_Pi_3D('symbolic_int')[
                                    np.ix_([idx],[V.index_of_weight(chi) for chi in tau.positive_weights(V)[y]], 
                                    [alpha.index_in_all_of_U(V.G) for alpha in gr_inv[y]])].sum(axis=0)
                                nrow,ncol=Mn.shape
                                #M=matrix(V.QV,Mn)
                                #print('type Ared',type(Ared))
                                #L0[0,col] = sum(
                                #    Mn[r,s] * Aredj.matrix_from_rows_and_columns(
                                #        [u for u in range(n_current) if u != r],
                                #        [v for v in range(n_current) if v != s]
                                #    ).det() * (-1)**(r + s)
                                #    for r in range(nrow)
                                #    for s in range(ncol)
                                #    if Mn[r,s] !=0
                                #)
                                # Alternative computation using Cramer
                                L0[0,col] = sum(
                                    Mn[r,s] * com_Aredj[r,s]
                                    for r in range(nrow)
                                    for s in range(ncol)
                                    if Mn[r,s] !=0
                                )                                
                        else : # In this case we use a symbolic method to avoid computation of all small minors.
                            print('A symbolic case with mult',mult_min,' and block of degree ',len(gr_inv[y]),'for:')
                            print(ineq)
                            Mn = V.T_Pi_3D('symbolic')[
                                    np.ix_([V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)],
                                    [V.index_of_weight(chi) for chi in tau.positive_weights(V)[y]], 
                                    [alpha.index_in_all_of_U(V.G) for alpha in gr_inv[y]])].sum(axis=0)
                            #print('matrix',Mn)
                            Jj=matrix(V.QV,Mn).det()
                            #print('Jj',Jj)
                            partial_derivatives: list[Polynomial] = [
                                Jj.derivative(V.QV.variable(chi))                                             
                                for chi in tau.orthogonal_weights(V)
                            ]
                            #print('partial der',partial_derivatives)
                            from sage.all import gcd as sage_gcd # type: ignore
                            #print('gcd',sage_gcd([Jj] + partial_derivatives))
                            Jjred: Polynomial = Jj//sage_gcd([Jj] + partial_derivatives)
                            Gradiant=[V.QV(Jjred.derivative(V.QV.variable(chi))) for chi in tau.orthogonal_weights(V)]
                            L0QV=matrix(V.QV,1,len(tau.orthogonal_weights(V)),Gradiant)
                            phi = V.T_Pi_3D(method_R0, 'dict')[p]
                            L0 = matrix(V.QZ, L0QV.apply_map(phi))
        
                        if L0*B0z_red*noyau != 0 :
                            #print('Bizarre')
                            return False 
    #print('et là ?')                    
    return True                    
                