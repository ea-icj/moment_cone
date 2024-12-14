from .rings import Matrix, vector, Polynomial

def Bezout_Inverse(LP):
    " LP is a list of polynomial pairwise coprime. Return a list of polynomial M_i such that Bezout inverse is sum(M_i*B_i)"
    res=[]
    for i,P in enumerate(LP):
        Q=prod(LP[:i]+LP[i+1:])
        u,A,B=xgcd(Q,P)
        res.append(A*Q)
    return(res)

def Kernel_modulo_P(ring,M : Matrix,LP,LIB,d)->Polynomial: # LP and LIB (list inverse Bezou) list of polynomials
    """M is a square matrix with entries in K[X]. We know the its Kernel in K[X]/prod(LP) is a submodule of rank 1 of (K[X]/prod(LP))^n. 
       The function returns a generator of this kernel by using the Chinese Theorem."""
    
    # Step 1: reduce M modulo each element of LP and then compute the kernel
    Vreduced=[]
    Qring=[]
    for P in LP:
        Qring.append(ring.quotient(P))
        Mreduced=M.apply_map(lambda entry: Qring[-1](entry))
        Kernel=Mreduced.right_kernel().basis()
        if len(b)!=1:
            return "Erreur de rang"
        else:
            Vreduced.append(b[0])
    # We now apply Bezout inverse to each coordinate of Vreduced      
    v=vector(ring,len(LP))
    for i in range(M.nrows()):
        v[i]=sum(Vreduced[k][i].lift()*LIB[k] for k in range(len(LP)))
    return(v)    
