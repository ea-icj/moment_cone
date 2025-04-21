__all__ = (
    "dim_gen_stab_of_K",
)

from random import randint
import itertools
import math
import numpy as np
from numpy.typing import NDArray


from .typing import *
from .linear_group import *
from .weight import *
from .root import *
from .representation import *
from .tau import *
from .rings import matrix, Matrix, vector, Vector, QQ, ZZ, I, real_part, imag_part

# Reference implementation
def dim_gen_stab_of_K_vref(
        T: NDArray[Any],
        ListK: Optional[Sequence[int]] = None, 
        ListChi: Optional[Sequence[int]] = None
    ) -> int:
    """
    Recursive function associating an integer to a list of matrices.

    Arguments:
    - T : 3 dimensional numpy.tab containing the action of a basis of Lie(K) on V (viewed as a real real vectorspace).
    - these matrices are the images of a basis of Lie(K) in End(V)
    - ListK is a list of elements in the basis of Lie(K) encoding a subalgebra
    - ListChi is a list of elements in the basis of V encoding a subspace
    Returns: an integer.
    """
    #for j in ListChi:
    #    print(
    # Default values of ListK and ListChi correspond to the entire array T
    assert T.shape[1] == T.shape[2]

    dk,dV,dV = T.shape
    if ListK is None:
        ListK=list(range(dk))
    if ListChi is None:
        ListChi=list(range(dV))

    # Check if all matrices are zero that is V is the trivial representation
    dk = len(ListK)
    if all(T[k,i,j]==0 for k,i,j in itertools.product(ListK,ListChi,ListChi)):
        return dk

    
    n: int = len(ListChi)  # Size of the square matrices that is dimension of the represention as a real vector space
    # Create the vector v in the representation
    #v = vector(ZZ, [randint(-9,9) for i in range(n)])
    v = np.random.randint(-9, 10, size=n) # higher bound in excluded in Numpy

    # Construct the matrix M
    #M = matrix(QQ, n, dk, lambda i, k: sum([T[ListK[k],ListChi[i],ListChi[j]] * v[j] for j in range(n)]))
    M = matrix(QQ, (T[np.ix_(ListK, ListChi, ListChi)] * v).sum(axis=-1).T)

    #from sympy import Matrix as Matrix_sympy
    #Ms = Matrix_sympy(n, dk, lambda i, k: sum(T[ListK[k]][ListChi[i]][ListChi[j]] * v[j] for j in range(n)))
    #Bs_tmp, pivots = Ms.T.rref()

    # Echelon form of M.transpose() to computation modulo the image F of M
    B_tmp = M.transpose().echelon_form().rref() # reduced echelon form
    B: Matrix = B_tmp.matrix_from_rows(B_tmp.pivot_rows()) # Suppress zero rows
    List_Pivots: list[int] = B.pivots()

    # Dimension of V/F
    dQ = n-len(List_Pivots)

    # If V/F is trivial then we can conclude
    if dQ == 0 :
        return(dk-n)
    
    if len(List_Pivots)==0:
        return dim_gen_stab_of_K(T)
    # The images of the elements of the canonical bases indexed by i not in List_Pivots form a basis Bc of V/F 
    List_Not_Pivots=[i for i in range(B.ncols()) if i not in List_Pivots] # TODO : Complexité quadratique. On peut aller plus vite comme dans complement_of_coordinate avec first...    # The image of e_i for i in List_Pivots is the ith columns of N in the basis Bc
    N: Matrix = -B.matrix_from_columns(List_Not_Pivots).transpose()
    #print(B,B_tmp,B,N)
    # Compute the basis of the left kernel of M. That a bases of the stabilizer of v.
    kernel_basis: list[Vector] = M.right_kernel().basis()

    #kernel_basis = Ms.nullspace()

    dk_stab=len(kernel_basis)

    # Determine the set I to form a basis of V/F ## 
    
    List_B=[matrix(QQ,dQ,dQ) for i in range(dk_stab)]
    TK = np.empty((dk_stab,n,len(List_Not_Pivots)),dtype=object) 
    T_stab = np.empty((dk_stab,dQ,dQ),dtype=object)
    for k, L in enumerate(kernel_basis) :
        for j,nj in enumerate(List_Not_Pivots):
            for i,ni in enumerate(ListChi):
                #action of the k-th element of (K^{tau})^v: V'->V^{tau} where V' is a supplementary supspace of V^{tau}/K^{tau}v
                TK[k,i,j] = sum([L[s]*T[ListK[s],ni,ListChi[nj]] for s in range(dk)]) 
            #Split TK
            #print(List_Pivots,type(List_Pivots),len(List_Pivots))
            #nv_pivots=vector(QQ,len(List_Pivots), lambda i: QQ(TK[k,List_Pivots[i],j]))
            vvv=vector(QQ,n, [TK[k,i,j] for i in range(n)])
            nv_pivots=vector(QQ,len(List_Pivots))
            for i,ip in enumerate(List_Pivots):
                nv_pivots[i]=TK[k,ip,j]
            nv_quot=vector(QQ,dQ) #, lambda i: QQ(TK[k,List_Not_Pivots[i],j]))
            for i,inp in enumerate(List_Not_Pivots):
                nv_quot[i]=TK[k,inp,j]
            nv_quot+=N*nv_pivots
            for i in range(dQ):
                T_stab[k,i,j]=nv_quot[i]
            
    # Recursive call with the new list of matrices
    return dim_gen_stab_of_K_vref(T_stab)


# Optimized implementation
def dim_gen_stab_of_K(
        T: NDArray[Any],
        ListK: Optional[Iterable[int]] = None,
        ListChi: Optional[Iterable[int]] = None
    ) -> int:
    """
    Recursive function associating an integer to a list of matrices.

    Arguments:
    - T : 3 dimensional numpy.array containing the action of a basis of Lie(K) on V (viewed as a real real vectorspace).
    - these matrices are the images of a basis of Lie(K) in End(V)
    - ListK is a list of elements in the basis of Lie(K) encoding a subalgebra
    - ListChi is a list of elements in the basis of V encoding a subspace
    Returns: an integer.
    """
    assert T.shape[1] == T.shape[2]
    dk, dV, dV = T.shape

    if ListK is None:
        ListK = np.arange(dk)
    if ListChi is None:
        ListChi = np.arange(dV)
    ListK = np.asarray(ListK)
    ListChi = np.asarray(ListChi)

    # Check if all matrices are zero that is V is the trivial representation
    dk = len(ListK)
    if np.all(T[np.ix_(ListK, ListChi, ListChi)] == 0):
        return dk
    
    # Size of the square matrices that is dimension of the representation as a real vector space
    n: int = len(ListChi)

    # Create the vector v in the representation
    v = np.random.randint(-9, 10, size=n) # higher bound in excluded in Numpy

    # Construct the matrix M
    # M_{i, k} = \sum_{j=0}^n T_{ListK_k, ListChi_i, ListChi_j} * v_j
    # with n = len(ListChi)
    # Since indices order change (k,i -> i,k), we need to use the transposition:
    M = matrix(QQ, (T[np.ix_(ListK, ListChi, ListChi)] * v).sum(axis=-1).T)

    # Echelon form of M.transpose() to computation modulo the image F of M
    B_tmp = M.transpose().echelon_form().rref() # reduced echelon form
    B: Matrix = B_tmp.matrix_from_rows(B_tmp.pivot_rows()) # Suppress zero rows
    List_Pivots = np.asarray(B.pivots())
        

    # Dimension of V/F
    dQ = n - len(List_Pivots)

    # If V/F is trivial then we can conclude
    if dQ == 0 :
        return(dk-n)
    
    #in the case where no Pivots where found, the random element lies in V^K (of dimension <dim V), and the slice theorem makes doesn't help in this case. So we restart the computation
    if len(List_Pivots)==0:
        return dim_gen_stab_of_K(T) 

    
    # The images of the elements of the canonical bases indexed by i not in List_Pivots form a basis Bc of V/F 
    List_Not_Pivots = np.array([i for i in range(B.ncols()) if i not in List_Pivots])

    # The image of e_i for i in List_Pivots is the ith columns of N in the basis Bc
    N: Matrix = -B.matrix_from_columns(List_Not_Pivots).transpose()

    # Converting N to an integer valued matrix
    N_lcm = math.lcm(*(int(Nij.denom()) for Ni in N for Nij in Ni if Nij != 1))
    N_int = np.empty((N.nrows(), N.ncols()), dtype=object) # dtype=object to use arbitrary precision integer of Python
    for i in range(N.nrows()):
        N_int[i, :] = [int(N_lcm * N[i, j]) for j in range(N.ncols())]

    # Compute the basis of the left kernel of M. That a bases of the stabilizer of v.
    kernel_basis: list[Vector] = M.right_kernel().basis()

    # Converting the basis to integer vectors
    kernel_basis_int = np.empty((len(kernel_basis), dk), dtype=object) # dtype=object to use arbitrary precision integer of Python
    for i, v in enumerate(kernel_basis):
        denom = math.lcm(*(frac.denom() for frac in v))
        kernel_basis_int[i, :] = [int(vi * denom) for vi in v]

    # Determine the set I to form a basis of V/F

    # Action of the k-th element of (K^{tau})^v: V'->V^{tau} where V' is a supplementary supspace of V^{tau}/K^{tau}v
    #TK[k,i,j] = sum([L[s]*T[ListK[s],ni,ListChi[nj]] for s in range(dk)]) 
    # TK_{k, i, j} = \sum_{l=0}^dk T_{ListK_l, ListChi_i, ListChi_{ListNotPivots_j}} * KB_{k, l}
    # with dk = len(ListK) and KB = kernel_basis_int.
    # Adding None placeholder to match indices order (k,l,i,j) in the sum, we get:
    # TK_{k, i, j} = \sum_{l=0}^dk T_{None, ListK_l, ListChi_i, ListChi_{ListNotPivots_j}} * KB_{k, l, None, None}
    # Thus:
    TK = (T[np.ix_(ListK, ListChi, ListChi[List_Not_Pivots])][None, :, :, :] * kernel_basis_int[:, :, None, None]).sum(axis=1)
    assert TK.dtype == object # ensuring that we still use Python int

    # nv_pivots_{k, i, j} = TK_{k, ListPivots_i, j}
    nv_pivots = TK[:, List_Pivots, :]

    # nv_quot_{k, i, j} = N_lcm TK_{k, ListNotPivots_i, j} + N_int @ nv_pivots_{k, *, j}
    # where * indicates the dimension of matrix@vector product and N_lcm is here to counterpart the integer conversion above.
    # In case a n-dimensional arrays (n > 2), numpy.matmul consider the first n-2 dimensions as batches
    # and apply a matrix product on the 2 trailing dimensions.
    nv_quot = N_lcm * TK[:, List_Not_Pivots, :] + np.matmul(N_int, nv_pivots)
    
    # T_stab is actually nv_quot...
    T_stab = nv_quot

    return dim_gen_stab_of_K(T_stab)

