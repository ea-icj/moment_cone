from .typing import *
from .dimension import *
from .weight import *
from .root import *
from .tau import *
from .rings import matrix, Matrix, vector, Vector, QQ, I, real_part, imag_part
from random import randint



def mat_C_to_R(M : Matrix) -> Matrix :
    "M is a matrix with complex coefficients. Replace each coefficient coefficien a+bI by a 2x2-matrix [a,-b,b,a]"
    A = real_part(M)
    B = imag_part(M)
    p = M.nrows()
    q = M.ncols() 
    R = matrix(QQ,2*p,2*q)
    for i in range(p):
        for j in range(q):
            R[2*i,2*j]=A[i,j]
            R[2*i+1,2*j+1]=A[i,j]
            R[2*i+1,2*j]=B[i,j]
            R[2*i,2*j+1]=-B[i,j]
    return(R)

def Lie_action_as_matrices_V(d : Dimension) -> list[Matrix]: # Remplace t-on par une méthode de d ?
    D=d.dimV 
    BaseK=Root.all_of_K(d)
    L=[]
    for beta in BaseK :
        k,i,j = beta.k,beta.i,beta.j
        # Matrix of the element [k,i,j] of the bases of Lie(K)
        M = matrix(QQ[I],D,D)
        d1 = Dimension(d[:k] + d[k+1:])
        for w in Weight.all(d1) : #list(itertools.product(*(range(di) for di in d1))):
            wj = list(w[:k])+[j]+list(w[k:])
            idj=Weight(wj).index_in(d)
            if i==j:
                M[idj,idj]=I
            else :    
                wi = list(w[:k])+[i]+list(w[k:]) # Eij. wj = wi
                idi=Weight(wi).index_in(d)
                if i<j :
                    M[idi,idj]=1
                    M[idj,idi]=-1
                else :
                    M[idi,idj]=I
                    M[idj,idi]=I
        L.append(M)    
    return(L)

def Lie_action_as_matrices_Vtau(tau : Tau,matrices) -> list[Matrix]: # matrices is a list of matrices
    d=tau.d
    Indices_V_tau=[chi.index_in(d) for chi in tau.orthogonal_weights]
    n=len(Indices_V_tau)
    res=[]

    # Todo faire tau.all_orthogonal_roots pour éviter les 3 cas ci-dessous
    for beta in tau.orthogonal_all_roots:
        res.append(matrices[beta.index_in_all_of_K(d)].matrix_from_rows_and_columns(Indices_V_tau,Indices_V_tau))
    return(res)


# Rename, like dim_of_stabilizer_in_K
def dim_gen_stab_of_K(matrices) -> int:
    """
    Recursive function associating an integer to a list of matrices.

    Arguments:
    - matrices: list [A_0, ..., A_{k-1}] of square matrices of size nxn over Q[i]
    - these matrices are the images of a basis of Lie(K) in End(V)
   
    Returns: an integer.
    """
    # Check if all matrices are zero that is V is the trivial representation
    dk=len(matrices) # dimension of K
    n = matrices[0].nrows()
    if all(A.is_zero() for A in matrices):
        #print('trivial representation of dimension',n,'dim K',dk)
        return dk

    n = matrices[0].nrows()  # Size of the square matrices
    # Create the vector v in the representation
    v = vector(QQ, [randint(-3,3) for i in range(n)])
    # Construct the matrix M
    M = matrix(QQ, n, dk, lambda i, k: sum([matrices[k][i,j] * v[j] for j in range(n)]))
    
    if M.is_zero(): # If M==0, we change v one for which M is not sero. Almost never used.
        for A in matrices:
            if not A.is_zero():
                for k,C in enumerate(A.columns()):
                    if not C.is_zero():
                        break
                break
        v = vector(QQ,n)
        v[k]=1
        M = matrix(QQ, n, dk, lambda i, k: sum([matrices[k][i,j] * v[j] for j in range(n)]))

    # Echelon form of M.transpose() to computation modulo the image F of M
    B_tmp = M.transpose().echelon_form().rref() # reduced echelon form
    B: Matrix = B_tmp.matrix_from_rows(B_tmp.pivot_rows()) # Suppress zero rows
    List_Pivots: list[int] = B.pivots()

    # Dimension of V/F
    qn = n-len(List_Pivots)

    # If V/F is trivial then we can conclude
    if qn == 0 :
        return(dk-n)
    
    # The images of the elements of the canonical bases indexed by i not in List_Pivots form a basis Bc of V/F 
    List_Not_Pivots=[i for i in range(B.ncols()) if i not in List_Pivots] # TODO : Complexité quadratique. On peut aller plus vite comme dans complement_of_coordinate avec first...    # The image of e_i for i in List_Pivots is the ith columns of N in the basis Bc
    N: Matrix = -B.matrix_from_columns(List_Not_Pivots).transpose()
  
    # Compute the basis of the left kernel of M. That a bases of the stabilizer of v.
    kernel_basis: list[Vector] = M.right_kernel().basis()
    dk_stab=len(kernel_basis)
    
    # Determine the set I to form a basis of V/F ## 
    
    List_B=[matrix(QQ,qn,qn) for i in range(dk_stab)] 
    for k, L in enumerate(kernel_basis) :
        for j,nj in enumerate(List_Not_Pivots):
            nv = cast(list[Vector], sum([L[i]*matrices[i].column(nj) for i in range(dk)]))
            
            #Split nv
            nv_pivots=vector(QQ,len(List_Pivots))
            nv_quot=vector(QQ,qn)
            ip=0;iq=0
            for i in range(n):
                if i in List_Pivots:
                   nv_pivots[ip]=nv[i]
                   ip+=1
                else :
                    nv_quot[iq]=nv[i]
                    iq+=1
                   
            nv_quot+=N*nv_pivots        
            for i in range(qn):
                List_B[k][i,j]=nv_quot[i]
            
    # Recursive call with the new list of matrices
    return dim_gen_stab_of_K(List_B)
