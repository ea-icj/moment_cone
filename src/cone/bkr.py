import numpy as np
from numpy.typing import NDArray
import itertools

from sage.all import SymmetricFunctions # type: ignore

from .typing import *
from .dimension import Dimension
from .tau import Tau
from .rings import ZZ, QQ, Vector, matrix

_sym_f = SymmetricFunctions(QQ).s()

def Kron_multi(L) -> int :
    """
    L is a list of partitions of the same weight, length at least 2
    return the multiple Kronecker coefficient
    """
    
    if len(L)==2:
        if L[1]==L[0]:
            return 1
        else:
            return 0
        
    product = _sym_f(L[0]).kronecker_product(_sym_f(L[1]))
    #### This part is unuseful but should go slightly faster with
    if len(L)==3 : # We look for L[2] in product
        for monomial, coeff in product.monomial_coefficients().items():
            if monomial==L[2]:
                return(coeff)
        return(0)    
    #### end unuseful
    
    # recursive part
    tot=0
    for monomial, coeff in product.monomial_coefficients().items():
        tot+=coeff*Kron_multi(L[2:]+[monomial])
    return(tot)    
        

def LR_multi(L,nu):  
    """
    L is a list of partitions, nu is a partition
    return the multiple LR-coeffient
    """
    import sage.libs.lrcalc.lrcalc as lr # type: ignore

    if len(L)==0:
        if nu==[]:
            return 1
        else :
            return 0
    if len(L)==1:
        if list(nu)==L[0]:
            return 1
        else:
            return 0
    if len(L)==2:
        
        return lr.lrcoef(nu,L[0],L[1])
    dict_mu=lr.mult(L[-2],L[-1])
    tot=0
    for mu in dict_mu.keys():
        tot+=dict_mu[mu]*LR_multi(L[:-2]+[mu],nu)
    return tot


def Kron_mat(Lambda):
    """
    Lambda is a table of partitions.
    return the product of the multiple Kronecker coeffients associated to the rows
    """
    p,s = Lambda.shape                  
    tot=1
    for i in range(p) :
        k=Kron_multi(Lambda[i,:])
        if k==0 :
            return(0)
        else :
            tot*=k
    return(tot)

def fct_weights_of_Nu(Nu): # Nu is a partial matrix with Partitions as entries
    """ 
        Nu is a partial matrix (as a list of columns) with Partitions as entries.
        This function replace each partition p by its weight, that is sum(p)
    """
    p,s = Nu.shape
    M = matrix(ZZ,p,s)
    for k in range(s):
        i=0
        while i < p and Nu[i,k] is not None:
            M[i,k]=sum(Nu[i,k])
            i+=1
    return M


def chi2Numat(chi, mult_tau ): # chi is a sage vector - mult_tau est un tau.reduced.mult
    """ Constructing nu by gluing consecutive elements of chi depending on given multiplicity.
    The result is a partial matrix with partitions (this is a property of chi due to the mathematical construction) as entries.
    Nu is a list of columns that are lists of partitions. Each partition is a list itself. 

    """
    from .utils import trim_zeros

    Nu = np.empty((max([len(x) for x in mult_tau]), len(mult_tau)), dtype=object)
    shift = 1
    for k,col in enumerate(mult_tau.blocks):
        for i,mult in enumerate(col):
            p=chi[shift:shift + mult]
            # Add in nu[-1] after removing tailing 0's
            Nu[i,k]=trim_zeros(p)
            shift += mult

    assert shift == len(chi)
    return Nu

def ineqs_in_Nl(l:int) : # the output is a list of inequalities for a Sagemath polytope in Z^l
  """
  Create the list of inequalities (tuples) saying that the polyhedron is contained in N^l. Namely
  0 1 0 0 0 ...
  0 0 1 0 0 ...
  0 0 0 1 0 ...
  0 0 0 0 1 ...
  ...
  """
  ineqs=[]
  v=[0]*(l+1)
  for i in range(l):
    v[i+1]=1
    ineqs.append(tuple(v))
    v[i+1]=0
  return ineqs



def Enumerate_delta(delta : int,ListP,sizenu) : # output : of list of integer vectors. sizenu a partial matrix of int as list of columns. ListP ??
    """ 
    return the list of possible weight vectors delta. That is satisfying condition 2 of Proposition prop:multCVtau.
    delta is a vector of nonnegative integers of size len(ListP).
    """
    # Polyhedron in Z^l
    l=len(ListP)
    # Polyhedron in N^l
    ineqs=ineqs_in_Nl(l)
    # Polyhedron in the hyperplan sum delta(j)=delta
    v=[-delta]+[1]*l
    eqs=[tuple(v)]

    # equalities corresponding to conditions 2
    for k in range(sizenu.ncols()): # Run over columns of Nu
        for i in range(sizenu.nrows()): # Run over the rows of the column
            if sizenu[i,k]>0:
                v=[-sizenu[i,k]]
                for j in range(l):
                    if ListP[j][k]==i:
                        v.append(1)
                    else:
                        v.append(0)
                eqs.append(tuple(v))

    #Create the polyhedron
    from sage.all import Polyhedron # type: ignore     
    PP=Polyhedron(ieqs=ineqs,eqns=eqs)

    #Return the integral points
    return [PP.dim(),PP.integral_points()]

def Fill_Table_of_Lambdas(delta, ListP, tau: Tau, d: Dimension, table: NDArray[Any]): # delta is a sage vector in ZZ -  output : void
    """
    Fill the table with the possible entries for Lambda.
    """
    # FIXME: why not reusing Partition from our library?
    from sage.all import Partitions # type: ignore

    for k in range(len(d)):
        for j in range(len(delta)):
            # max_length triggers the type checker.
            table[j,k]=Partitions(delta[j], max_length=tau.components[k][ListP[j][k]]).list() # type: ignore
    

def LR_mat(Nu,Lambda,ListP) -> int : # Lambda is a table of partitions
                       # Nu is a "partial table" of partitions
                       # ListP is a list of list of row indexes
    # return the product of the multiple LR-coeffients associated to the rows

    tot=1
    p,s = Nu.shape
    for k in range(s):
        i=0
        while i < p and Nu[i,k] is not None: 
            Lambda_m=[]  
            for j in range(len(ListP)):
                if ListP[j][k]==i:
                    Lambda_m.append(Lambda[j][k])
          
            lrm=LR_multi(Lambda_m,Nu[i,k])
            
            if lrm==0:            
                return(0)
            else:            
                tot*=lrm
            i+=1
            
    return(tot)        
                       

def Multiplicity_SV_tau(tau : Tau,chi : Vector) -> int:

    d=tau.d
    delta=chi[0]
    ListP=tau.indices_in_tau_red_that_sum_to_zero
    Nu=chi2Numat(chi,tau.reduced.mult)  # Nu is a dominant weight for G^\tau

    w_Nu=fct_weights_of_Nu(Nu)
    mult=0

    for delta_j in Enumerate_delta(delta,ListP,w_Nu)[1]: # Run over the dela_j satisfying Condition 2
        table_Lambda = np.empty((len(ListP), len(d)), dtype=object)
        Fill_Table_of_Lambdas(delta_j,ListP,Tau(tau.reduced.mult),d,table_Lambda) # The list of tabulars of partitions satisfying Condition 1 and 3, for the given delta_j

        # Create the list of Lambda's from table_Lambda
        
    
        #List_of_Lambdas = [np.array(p, dtype=object).reshape(table_Lambda.shape) for p in itertools.product(*table_Lambda.ravel())]
        # reshape ne marche pas si p est une liste de listes à 1 éléments. Donc à la main:
        
        List_of_Lambdas=[]
        a,b=table_Lambda.shape
    
        for p in itertools.product(*table_Lambda.ravel()):
            T= np.empty(table_Lambda.shape, dtype=object)
            for i in range(a):
                for j in range(b):
                    T[i,j]=p[i*b+j]
            List_of_Lambdas.append(T)
            
        
            
        for Lambda in List_of_Lambdas : 
            
            LR_coeff=LR_mat(Nu,Lambda,ListP)
           
            if LR_coeff!=0:
                K_coeff=Kron_mat(Lambda)
                mult+=LR_coeff*K_coeff
    return(mult)            

def Is_Multiplicity_SV_tau_one(tau : Tau,chi : Vector) -> bool:

    d=tau.d
    delta=chi[0]
    ListP=tau.indices_in_tau_red_that_sum_to_zero
    Nu=chi2Numat(chi,tau.reduced.mult)  # Nu is a dominant weight for G^\tau

    w_Nu=fct_weights_of_Nu(Nu)
    mult=0
    Integral_Points=Enumerate_delta(delta,ListP,w_Nu)
    # Regular case: we check that the dimension of P is 0
    if tau.is_dom_reg :
        return(Integral_Points[0]==0)

    #Non regular case
    for delta_j in Integral_Points[1]: # Run over the dela_j satisfying Condition 2
        table_Lambda = np.empty((len(ListP), len(d)), dtype=object)
        Fill_Table_of_Lambdas(delta_j,ListP,Tau(tau.reduced.mult),d,table_Lambda) # The list of tabulars of partitions satisfying Condition 1 and 3, for the given delta_j

        # Create the list of Lambda's from table_Lambda
        
    
        #List_of_Lambdas = [np.array(p, dtype=object).reshape(table_Lambda.shape) for p in itertools.product(*table_Lambda.ravel())]
        # reshape ne marche pas si p est une liste de listes à 1 éléments. Donc à la main:
        
        List_of_Lambdas=[]
        a,b=table_Lambda.shape
    
        for p in itertools.product(*table_Lambda.ravel()):
            T= np.empty(table_Lambda.shape, dtype=object)
            for i in range(a):
                for j in range(b):
                    T[i,j]=p[i*b+j]
            List_of_Lambdas.append(T)
            
        
            
        for Lambda in List_of_Lambdas : 
            
            LR_coeff=LR_mat(Nu,Lambda,ListP)
           
            if LR_coeff!=0:
                K_coeff=Kron_mat(Lambda)
                mult+=LR_coeff*K_coeff
            if mult>1: return(False)
    return(mult==1)            
