from .rings import Matrix, vector, Polynomial

def Bezout_Inverse(LP,ring): #TODO : typer
    """ 
    LP is a list of polynomial pairwise coprime. 
    Return a list of polynomial M_i such that Bezout inverse is sum(M_i*B_i)
    """
    res=[]
    for i,P in enumerate(LP):
        Q=ring(1)
        for pol in LP[:i]+LP[i+1:]:
            Q*=pol
        u,A,B=Q.xgcd(P)
        res.append(A*Q)
    return(res)

def Kernel_modulo_P(ring,M : Matrix,LP,LIB)->Polynomial: # LP and LIB (list inverse Bezou) list of polynomials. #TODO : typer
    """
    M is a square matrix with entries in K[X]. 
    We know that the its Kernel in K[X]/prod(LP) is a submodule of rank 1 of (K[X]/prod(LP))^n. 
    The function returns a generator of this kernel by using the Chinese Theorem.
    """
    
    # Step 1: reduce M modulo each element of LP and then compute the kernel
    Vreduced=[]
    Qring=[]
    for P in LP:
        Qring.append(ring.quotient(P))
        Mreduced=M.apply_map(lambda entry: Qring[-1](entry))
        Kernel_basis=Mreduced.right_kernel().basis()
        if len(Kernel_basis)!=1:
            return "Erreur de rang"
        else:
            Vreduced.append(Kernel_basis[0])
    
    # We now apply Bezout inverse to each coordinate of Vreduced      
    v=vector(ring,M.nrows())
    for i in range(M.nrows()):
        v[i]=sum(Vreduced[k][i].lift()*LIB[k] for k in range(len(LP)))
    return(v)    
