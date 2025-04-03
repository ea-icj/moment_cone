__all__ = (
    "dim_gen_stab_of_K",
)

from random import randint
import itertools
import numpy as np
import numpy.typing as npt


from .typing import *
from .linear_group import *
from .weight import *
from .root import *
from .representation import *
from .tau import *
from .rings import matrix, Matrix, vector, Vector, QQ, ZZ, I, real_part, imag_part





def dim_gen_stab_of_K(T,ListK = None,ListChi = None) -> int: # New
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
    dk,dV,dV = T.shape
    if ListK == None:
        ListK=list(range(dk))
    if ListChi == None:
        ListChi=list(range(dV))

    # Check if all matrices are zero that is V is the trivial representation
    dk = len(ListK)
    if all(T[k,i,j]==0 for k,i,j in itertools.product(ListK,ListChi,ListChi)):
        return dk

    
    n: int = len(ListChi)  # Size of the square matrices that is dimension of the represention as a real vector space
    # Create the vector v in the representation
    v = vector(ZZ, [randint(-9,9) for i in range(n)])
    # Construct the matrix M
    M = matrix(QQ, n, dk, lambda i, k: sum([T[ListK[k],ListChi[i],ListChi[j]] * v[j] for j in range(n)]))
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
    return dim_gen_stab_of_K(T_stab)

