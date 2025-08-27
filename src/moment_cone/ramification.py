__all__ = (
    'compute_d_minus_2_minors_via_adjugate',
    'is_not_contracted',
    'Is_Ram_contracted',
)

from random import randint
import itertools
import numpy as np
from numpy.typing import NDArray
from sage.all import NumberField # type: ignore

from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .representation import *
from .inequality import *
from .permutation import *
from .rings import QQ, I, matrix, vector, Matrix, Polynomial, PolynomialRing
from .utils import prod,fl_dic

    
from itertools import combinations

def compute_d_minus_2_minors_via_adjugate(A: Matrix) -> dict[tuple[tuple[int, ...], tuple[int, ...]], int]:
    d = A.nrows()
    minors_dict: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    
    # Générer toutes les sous-matrices B en supprimant une ligne i et une colonne j
    for i in range(1,d):
        for j in range(1,d):
            # Sous-matrice B = A sans la ligne i et la colonne j (taille (d-1) x (d-1))
            B = A.matrix_from_rows_and_columns(
                [k for k in range(d) if k != i],
                [l for l in range(d) if l != j]
            )
            
            # Calcul de l'adjugate de B
            adj_B = B.adjugate().transpose()
            
            # Parcourir les coefficients de adj(B)
            for k in range(i):
                for l in range(j):
                    # Trouver les indices originaux dans A
                    # Lignes supprimées : i et la k-ème ligne restante (dans B)
                    #rows_in_A = [i] + [m for m in range(d) if m != i][k]
                    # Colonnes supprimées : j et la l-ème colonne restante (dans B)
                    #cols_in_A = [j] + [n for n in range(d) if n != j][l]
                    
                    # Le mineur est adj_B[k, l], avec un signe dépendant des positions
                    #sign = (-1)**(sum(rows_in_A) + sum(cols_in_A))
                    #minor = sign * adj_B[k, l]
                    
                    # Stocker dans le dictionnaire (en utilisant frozenset pour l'ordre)
                    I = (k,i)
                    J = (l,j)
                    minors_dict[(I, J)] = (-1)**(k+l)*adj_B[k,l] 
    
    return minors_dict

def taylor_det_psi(k: int, A: Matrix, psi: Matrix, h_symbolic: bool = True) -> Polynomial:
    """
    Compute the Taylor expansion of det(A + ψ(h)) up to order k_max,
    as a polynomial in symbolic variables h_0, ..., h_{s-1}.

    Args:
        k_max (int): Maximum order of the expansion.
        A (matrix): Matrix ψ(a) (constant term).
        psi (list): 3D array where psi[k][i][j] is the coefficient of h_k in ψ(h)[i,j].
        h_symbolic (bool): If True, use symbolic h_k (default=True).

    Returns:
        Dict[int, Polynomial]: Dictionary mapping orders k to their polynomial terms in h.
    """
    assert k >= 2
    d = A.nrows()  # Deduce size from A
    term = 0
        
    dic_minors: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    B: Matrix = A
    if k == 2 :
        dic_minors = compute_d_minus_2_minors_via_adjugate(A)

    for rows in combinations(range(d), k):
        if k > 2 : 
            B = A.delete_rows(rows)
        for cols in combinations(range(d), k):    
            
            # Minor of A (delete rows/cols)
            if k > 2: 
                minor_A = B.delete_columns(cols).det()
            else :
                minor_A = dic_minors[(rows,cols)] 
        
            if minor_A != 0 :     
                # Product of H entries: ψ(h)[i1,j1] * ... * ψ(h)[ik,jk]
                h_monomial = psi[rows,cols].det() 
                
                # Sign: (-1)^(sum of rows + cols)
                sign = (-1)**(sum(rows) + sum(cols))
                term += sign * minor_A * h_monomial
    return term
    

  
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

    >>> from moment_cone import *
    >>> G = LinearGroup((2, 3, 4, 1))
    >>> V = KroneckerRepresentation(G)
    >>> tau = Tau.from_flatten([6, 2, 1, 4, 1, 4, 5, 3, 1, 1], G)
    >>> tau
    6 2 | 1 4 1 | 4 5 3 1 | 1
    >>> roots = list(Root.all_of_U(G))[:10]
    >>> positive_weights = list(itertools.chain.from_iterable(tau.positive_weights(V).values()))
    >>> non_positive_weights = list(itertools.chain.from_iterable(tau.non_positive_weights(V).values()))
    >>> is_not_contracted(roots, V, "probabilistic", positive_weights, non_positive_weights)
    False
    >>> is_not_contracted(roots, V, "symbolic", positive_weights, non_positive_weights)
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
    from .utils import getLogger
    logger = getLogger("ramification.Is_Ram_contracted")

    ws=ineq.w
    Inv_w=list(root for root in ineq.inversions)
    dU=len(Inv_w)
    if dU<=1 : 
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
    List_C_mult_dec=sorted(gr_inv.keys(), reverse=True) 
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
        for x in List_C_mult_dec : 
            n_block=len(gr_inv[x])
            Azi=Az.submatrix(shift,shift,n_block,n_block)
            
            Ji=Azi.det()
            Azi_com=Azi.adjugate()
            List_Comatices.append(Azi_com.transpose())
            Jz*=Ji
            Blocks_Az.append(Azi)
            Factors_Jz.append(Ji)
            List_deltas.append(dict(Ji.factor()))
            shift+=n_block
        
        # Run over the diagonal blocks of A
        List_checked_deltas=[]
        for i,x in enumerate(List_C_mult_dec) :
            # Current block and its determinant
            Azi=Blocks_Az[i]
            Ji=Factors_Jz[i]
            n_current=len(gr_inv[x])

            for delta, m in List_deltas[i].items():
                if delta not in List_checked_deltas:   
                    List_checked_deltas.append(delta)
                                          
                    
                    K=NumberField(delta, 'a')
                    a = K.gen()
                    Ared=Az.matrix_from_columns(range(sizeblocks[i], Az.ncols())).change_ring(K)
                    # Block used to compoute L0
                    
                    if Jz % delta**2 != 0  or Ared.rank() == Ared.ncols()-1 :
                        noyau=Ared.right_kernel().basis()[0] 
                        
                        ### Computation of L0
                        # choose a block of minimal multiplicity. We first list the occurences of delta
                        occ_delta=[[i,m,x]]
                        for j,y in enumerate(List_C_mult_dec[i+1:]):
                            if delta in List_deltas[i+j+1].keys() :
                                occ_delta.append([j+i+1,List_deltas[j+i+1][delta],y])
                        
                        jmin,mult_min,y_weight = min(occ_delta, key=lambda occ: [occ[1],len(gr_inv[occ[2]])])
                        
                        # We minimize multiplicity of delta. If several ones, we minimize the size of the block
                           
                        if mult_min==1:
                            nA = Ared.nrows()
                            L0=matrix(K,1,len(zw_idx))
                            Aredj=Blocks_Az[jmin].change_ring(K) # block jmin modulo delta
                            com_Aredj=List_Comatices[jmin].change_ring(K) # its comatrix
                            for col,idx in enumerate(zw_idx):
                                Mn = V.T_Pi_3D('symbolic_int')[
                                    np.ix_([idx],[V.index_of_weight(chi) for chi in tau.positive_weights(V)[y_weight]], 
                                    [alpha.index_in_all_of_U(V.G) for alpha in gr_inv[y_weight]])].sum(axis=0)
                                nrow,ncol=Mn.shape
                                L0[0,col] = sum(
                                    Mn[r,s] * com_Aredj[r,s]
                                    for r in range(nrow)
                                    for s in range(ncol)
                                    if Mn[r,s] !=0
                                )
                                                        
                        else : # In this case we use a symbolic method to avoid computation of all small minors.
                            logger.debug(f'A symbolic case with mult {mult_min} and block of degree {len(gr_inv[y_weight])} for:')
                            logger.debug(str(ineq))
                            
                            Mn = V.T_Pi_3D('symbolic')[
                                    np.ix_([V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)],
                                    [V.index_of_weight(chi) for chi in tau.positive_weights(V)[y_weight]], 
                                    [alpha.index_in_all_of_U(V.G) for alpha in gr_inv[y_weight]])].sum(axis=0)
                            
                            psi = matrix(V.QV, Mn)
                            
                            Aredj=Blocks_Az[jmin].change_ring(K)
                            
                            taylor_term=taylor_det_psi(mult_min, Aredj, psi)
                            
                            # Computation of L0 such that L0^mult_min = constant * taylor_term
                            variables = [V.QV.variable(chi) for chi in tau.orthogonal_weights(V)]
                            vars = taylor_term.variables()
                            var0 = vars[0]
                            coef_norm = taylor_term.coefficient({var0: mult_min})
                            taylor_term = taylor_term/coef_norm 
                            coeffs = [1]
                            for var in vars[1:]:
                                coeffs.append(taylor_term.coefficient({var0: mult_min-1,var:1})/mult_min)
                            coeffs_extended = [dict(zip(vars,coeffs)).get(var, 0) for var in variables]
                            L0 = vector(K, coeffs_extended)

                        B0z_red=B0z.matrix_from_columns(range(sizeblocks[i], Az.ncols())).change_ring(K)
                        
                        if (L0*B0z_red*noyau)[0] != 0 :
                            return False 
                  
    return True                    
                