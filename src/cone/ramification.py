from .typing import *
from .tau import Tau
from .root import Root
from .weight import Weight
from .cone_dimension import point_vect, vector, matrix
from .spaces import action_op_el
from .inequality import *
from .permutation import *
from .kx_mod import *
from random import randint
import itertools

# FIXME: we get d from tau but in the current code, it will leads to recreate the rings for each tau.
def is_not_contracted(inversions_v: Sequence[Root], tau: Tau, method: Method) -> bool:
    """
    ???

    The vector v is chosen depending on the method: probabilistic or symbolic.

    FIXME: example

    Example:
    >>> from cone import *
    >>> d = Dimension((2, 3, 4))
    >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
    >>> tau
    1 | 6 2 | 1 4 1 | 4 5 3 1
    >>> roots = list(Root.all_of_U(d))[:10]
    >>> is_not_contracted(roots, tau, "probabilistic")
    False
    >>> is_not_contracted(roots, tau, "symbolic")
    False
    """
    d = tau.d

    # Ring depending on the computational method
    if method == "probabilistic":
        ring = d.QI
    elif method == "symbolic":
        ring = d.QV
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    
    # FIXME: do we stay we list conversion at each call?
    # Maybe grading root and weight should be implemented using a more convenient class?
    from itertools import chain
    non_positive_weights = list(chain.from_iterable(tau.non_positive_weights.values())) # todo : Cela est aussi fait avant l'appel à la fonction.
    positive_weights = list(chain.from_iterable(tau.positive_weights.values())) # todo : Cela est aussi fait avant l'appel à la fonction.

    # FIXME: why not positive_weights instead of non_negative_weights?
    v = point_vect(non_positive_weights, d, ring, bounds=(-1000, 1000))
    list_inversions_v=list(inversions_v)
    A = matrix(ring, len(positive_weights), len(list_inversions_v))
    for j, root in enumerate(list_inversions_v):
        uv = action_op_el(root, v, d)
        for i, chi in enumerate(positive_weights):
            A[i, j] = uv[chi.index_in(d)]

    rank_A = A.change_ring(ring.fraction_field()).rank()
    return rank_A == len(list_inversions_v)

def Compute_JA(ineq : Inequality) : # Return a dictionnary polynom :  int
    tau=ineq.tau
    d = tau.d
    ring = d.QV
    # a generic vector in VV^tau
    zero_weights = tau.orthogonal_weights
    v = point_vect(zero_weights, d, ring, bounds=(-100, 100)) # bounds unuseful here
    # inversions of w
    Inv_w=ineq.inversions
    #gr = grading_dictionary(ineq.inversions, tau.dot_root)
    gr = tau.grading_roots_in(ineq.inversions)
    J={}
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau        
        M=matrix(ring,len(gr[x]))
        for col,root in enumerate(gr[x]): # List of roots such that tau.scalar(root)=x
            uv=action_op_el(root, v, d)
            for row, chi in enumerate(tau.positive_weights[x]): # List of weights such that tau.scalar(chi)=x 
                M[row,col]=uv[chi.index_in(d)]
        #print('M',M)        
        Jb=dict(M.det().factor())   
        for F in Jb.keys(): # We could make a function add_dictionaries
            if F in J.keys():
                J[F]+=Jb[F]
            else:
                J[F]=Jb[F]
    return(J)

def Smith_n_1(A):
    "Compute the gcd of the minors of A of size n-1"
    combinaisons = list(itertools.combinations(range(A.nrows()),A.nrows()-1))
    ring=A.base_ring()
    pgcd = ring(0)
    
    # Run over the n-1 x n-1 minors
    for lignes in combinaisons:
        for colonnes in combinaisons:
            # Extract submatrices
            sous_matrice = A.matrix_from_rows_and_columns(lignes, colonnes)
            
            # Det            
            det_Smith = sous_matrice.det()
            
            # Upgrade GCD
            pgcd = pgcd.gcd(det_Smith)
            
            #If the gcd is 1 stop
            if pgcd == 1:
                return pgcd    
    return pgcd

def Is_Ram_contracted(ineq : Inequality, method_S: Method, method_R0: Method) -> bool :
    d=ineq.tau.d
    if method_S == "probabilistic": 
        ring = d.QI  # ou QQ ?
        ring2= d.QZ
    elif method_S == "symbolic":
        ring = d.QV
        ring2= d.QV2
    else:
        raise ValueError(f"Invalid value {method} of the computation method")
    
    ws=ineq.w
    tau=ineq.tau
    Inv_w=list(root for root in ineq.inversions)
    dU=len(Inv_w)
    
    # Creation of sorted lists of weights todo: modifier en utilisant mieux la classe tau
    Neg0_Weights_dic=tau.non_positive_weights
    Neg0_Weights_sorted=[]
    for x in sorted(Neg0_Weights_dic.keys(),reverse=True):
        Neg0_Weights_sorted+=Neg0_Weights_dic[x]
    Pos_Weights_dic=tau.positive_weights
    Pos_Weights_sorted=[]
    for x in sorted(Pos_Weights_dic.keys(),reverse=True):
        Pos_Weights_sorted+=Pos_Weights_dic[x]
    if dU==0 : # Trivial case. Can we replace by <=1
        return(True)

    ### Divisors of the boudary
    for k,w in enumerate(ws):
        for v in w.covering_relations_strong_Bruhat:
            if v.is_min_rep(tau.reduced.mult[k]): 
                vs=Permutation(list(ws[:k])+[v]+list(ws[k+1:]))
                ineqv=Inequality(tau,vs)
                if is_not_contracted(ineqv.inversions,tau,method_S) :
                    return(False)

    #return(True) # Freeze the second checking
    ### Divisor R_0
    Jf=Compute_JA(ineq) # The Jacobian factorized as a dictionnary
    #print('Jf',Jf)
    #print(list(Jf.keys())[0] in d.QV)
    J_square_free=1
    for pol in Jf.keys():
        J_square_free*=pol # todo : prod(list(Jf.keys())) ne semble pas fonctionner
    if len(Jf.keys())!=len(dict(J_square_free.factor()).keys()):
        print('Error in factor with:',Jf,J_square_free)
    #print('J square free',J_square_free)        
    # gradiant of J_square_free
    L0=matrix(d.QV,1,len(Neg0_Weights_dic[0]))
    for col,chi in enumerate(Neg0_Weights_dic[0]) :
        L0[0,col]=J_square_free.derivative(d.QV.variable(chi))
    #L=matrix(d.QV,1,d.dimV)# A suprimer
    #for col,chi in enumerate(Weight.all(d)) : #A suprimer
    #    L[0,col]=J_square_free.derivative(d.QV.variable(chi))#A suprimer

    #print('L',L)
    #print('L0',L0)
    #print(len(L))
    #cpt=0
    #for i in range(27):
    #    if L[0,i]!=0 : cpt+=1
    #print('nb var in J',cpt)
    
    # Generic point v of V(tau<=0) and matrix of Tpi at (e,v)
    v = point_vect(Neg0_Weights_sorted,d,d.QV, bounds=(-4, 4)) # Augmenter les bornes
    #v = point_vect(Neg0_Weights_sorted,d,d.QV, bounds=(-1000, 1000))
    T=matrix(d.QV,d.dimV,dU)
    gr = tau.grading_roots_in(ineq.inversions)
    col=0
    for x in sorted(gr.keys(),reverse=True): # Choose a diagonal block of Tpi that is a weight of tau        
        for root in gr[x]: # List of roots such that tau.scalar(root)=x
            uv=action_op_el(root, v, d)
            row=0
            for row,chi in enumerate(Pos_Weights_sorted+Neg0_Weights_sorted):
                T[row,col]=uv[chi.index_in(d)]
            col+=1
    A=T[:dU,0:dU] # equivalent to T.submatrix(0, 0, shift, shift)
    B0=T[dU:dU+len(Neg0_Weights_dic[0]),0:dU] #0:dU is any
    
    # The line
    subs_dict = {}    
    for chi in Neg0_Weights_sorted:
        if method_R0 == "probabilistic":
            subs_dict[d.QV.variable(chi)]= randint(-500,500)*d.QZ('z')+randint(-500,500)# todo :Tester l'effet du changement de 500. Doit-on mettre du I ? 
        else:
            subs_dict[d.QV.variable(chi)]= d.QV2.variable(chi)*d.QV2('z')+d.QV2.variable(chi) # Mettre a et b 
    # Substitutions
    Az=A.subs(subs_dict)
    #D, U, V = Az.smith_form()
    #print('Smith')
    #for i in range(D.nrows()):
    #    print(D[i,i])
    B0z=B0.subs(subs_dict)
    L0z=L0.subs(subs_dict)
    Jz=J_square_free.subs(subs_dict)


    # Computation of reduced delta as factorized polynomial
    Smith_n_un=Smith_n_1(Az)
    s=Jz.degree()-Smith_n_un.degree()
    exponent=max(s-1,1)
    pgcd=Jz.gcd(Smith_n_un**exponent)
    #print('Jz,Smith_n_un,pgcd:',Jz,Smith_n_un,pgcd)
    delta1 = Jz // pgcd # delta1 =Jz.quo(pgcd)
    delta = delta1 // delta1.gcd(delta1.derivative()) #delta = delta1.quo(gcd(delta1,delta_1.derivative()))
    Ldelta=list(dict(delta.factor()).keys())
    #print('delta',delta,Ldelta)

    # Computation of Bezout inverse
    LIB=Bezout_Inverse(Ldelta,d.QZ)
    
    # Kernel of Az modulo delta
    noyau=Kernel_modulo_P(ring2,Az,Ldelta,LIB,d)

    #print('\n noyau\n',noyau,'\n\n')
    #for i in range(Az.nrows()):
    #    print((Az*noyau)[i] % delta)
    # Check divisibility
    check=L0z*B0z*noyau
    check=d.QZ(check[0])
    #print('check',check)
    quo, rem = check.quo_rem(delta)
    #print('A Quotient:',rem)
    return(rem==0)
   
