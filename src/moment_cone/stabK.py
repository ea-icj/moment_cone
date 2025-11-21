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


def pivot_columns_rref(M: Matrix) -> list[int]:
    """
    Returns the indices of the pivot columns (where the 1's are) in a scaled reduced-form matrix,
 assumed to have no zero rows and maximum rank.
    """
    nrows = M.nrows()
    pivots = []
    col = 0
    for row in range(nrows):
        while M[row, col] == 0:
            col += 1
        pivots.append(col)
        col += 1  # you can move forward directly because the pivots are strictly on the right line by line
    return pivots

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
    # Since indices order change (k,i -> i,k), we need to use the transposition:
    M = matrix(QQ, (T[np.ix_(ListK, ListChi, ListChi)] * v).sum(axis=-1).T)

    # Echelon form of M.transpose() to computation modulo the image F of M
    B_tmp = M.transpose().echelon_form().rref() # reduced echelon form
    #B: Matrix = B_tmp.matrix_from_rows(B_tmp.pivot_rows()) # Suppress zero rows
    rkB=B_tmp.rank(algorithm='flint')
    B = matrix(QQ, map(B_tmp.row, range(rkB)))
    List_Pivots = np.asarray(pivot_columns_rref(B))
        

    # Dimension of V/F
    dQ = n - rkB

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
    kernel_basis: list[Vector] = M.right_kernel(algorithm='flint').basis()

    # Converting the basis to integer vectors
    kernel_basis_int = np.empty((len(kernel_basis), dk), dtype=object) # dtype=object to use arbitrary precision integer of Python
    for i, v in enumerate(kernel_basis):
        denom = math.lcm(*(frac.denom() for frac in v))
        kernel_basis_int[i, :] = [int(vi * denom) for vi in v]

    # Determine the set I to form a basis of V/F

    TK = (T[np.ix_(ListK, ListChi, ListChi[List_Not_Pivots])][None, :, :, :] * kernel_basis_int[:, :, None, None]).sum(axis=1)
    assert TK.dtype == object # ensuring that we still use Python int

    # nv_pivots_{k, i, j} = TK_{k, ListPivots_i, j}
    nv_pivots = TK[:, List_Pivots, :]

    nv_quot = N_lcm * TK[:, List_Not_Pivots, :] + np.matmul(N_int, nv_pivots)
    
    # T_stab is actually nv_quot...
    T_stab = nv_quot

    return dim_gen_stab_of_K(T_stab)

