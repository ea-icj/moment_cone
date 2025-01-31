import numpy as np
import itertools

from sage.all import SymmetricFunctions,ZZ,QQ,vector,matrix,Polyhedron,Partitions
import sage.libs.lrcalc.lrcalc as lr

from .typing import *
from .utils import *
from .partition import *
from .group import *
from .weight import *
from .tau import *
from .rep import *
from .inequality import *

sym_f = SymmetricFunctions(QQ).s()


def Kron_multi_old(L) -> int :
    """
    L is a list of partitions of the same weight, length at least 2
    return the multiple Kronecker coeffient
    """
    
    if len(L)==2:
        if L[1]==L[0]:
            return 1
        else:
            return 0
        
    product = sym_f(list(L[0])).kronecker_product(sym_f(list(L[1])))
    #### This part is unuseful but should go slightly faster with
    if len(L)==3 : # We look for L[2] in product
        for monomial, coeff in product.monomial_coefficients().items():
            #print('monom',monomial,type(monomial),OurPartition(monomial))
            if OurPartition(list(monomial))==L[2]:
                return(coeff)
        return(0)    
    #### end unuseful
    
    # recursive part
    tot=0
    for monomial, coeff in product.monomial_coefficients().items():
        tot+=coeff*Kron_multi(L[2:]+[monomial])
    return(tot)    


def Kron_multi(L) -> int:
    """
    L is a list of partitions of the same weight, length at least 2
    return the multiple Kronecker coeffient
    """
    from sage.all import Partition as SagePartition
    product = sym_f(tuple(L[0]))
    for p in L[1:-1]:
        product = product.kronecker_product(sym_f(tuple(p)))
    return product.monomial_coefficients().get(SagePartition(L[-1]), 0)


def LR_multi(L,nu) -> int:
    
    """
    L is a list of partitions, nu is a partition
    return the multiple LR-coeffient
    """
    #print('LR_mult',L,nu,type(L))
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
        #print('L',L[:-2])
        #print(list(mu),type(mu))
        #print(L[:-2]+[list(mu)])
        tot+=dict_mu[mu]*LR_multi(L[:-2]+[mu],nu)
    return tot


def Kron_mat(Lambda) -> int:
    """
    Lambda is a rectangular table of partitions.
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

def LR_mat_fb(Mu,Nu) -> int:
    """
    fb means fermion of boson
    Mu is a rectangular table of partitions.
    Nu is a table of partitions with one column with the same number of rows as Lambda
    return the product of the multiple LR-coeffients associated to the rows
    """
    p,s = Mu.shape                  
    tot=1
    for j in range(s) :
        k=LR_multi([mu for mu in Mu[:,j]],Nu[j,0])
        if k==0 :
            return(0)
        else :
            tot*=k
    return(tot)

def Pleth_list(Lambda,Theta,Mu):
    """
    Lambda, Theta, Mu are lists of partitions
    Compute the product of the corresponding list of plethysms
    """
    res=1
    for l,t,n in zip(Lambda, Theta, Mu):
        pl=sym_f(list(l)).plethysm(sym_f(list(t))).coefficient(n)
        if pl == 0 :
            return 0
        else :
            res*=pl
    return res        
            
def fct_weights_of_Nu(Nu) -> matrix : # Nu is a partial matrix with Partitions as entries
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
    Nu = np.empty((max([len(x) for x in mult_tau]), len(mult_tau)), dtype=object)
    shift = 0
    for k,col in enumerate(mult_tau.blocks):
        for i,mult in enumerate(col):
            p=chi[shift:shift + mult]
            # Add in nu[-1] after removing tailing 0's
            Nu[i,k]=trim_zeros(p)
            shift += mult
            
    #print('chi',chi)
    #print('mult_tau',mult_tau)
    #print('Nu',Nu)
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



def Enumerate_delta(ListP : list[list[int]],sizenu,V : Representation,delta : int=0) : # output : of list of integer vectors. sizenu a partial matrix of int as list of columns. ListP ??
    """ 
    return the list of possible weight vectors delta. That is satisfying condition 2 of Proposition prop:multCVtau.
    delta is a vector of nonnegative integers of size len(ListP).
    """
    #print('sizenu',sizenu,sizenu.nrows())
    # Polyhedron in Z^l
    l=len(ListP)
    # Polyhedron in N^l
    ineqs=ineqs_in_Nl(l)
    if V.type == 'kron' : #TODO : modifier sans ccomponent
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
    else:
        eqs=[]
        for j in range(sizenu.nrows()):
        #for j,n in enumerate(sizenu[0]):
            eqs.append(tuple([-sizenu[j,0]]+[I[j] for I in ListP]))            
    #print('eqs',eqs)
    #print('ineq',ineqs)
    #Create the polyhedron        
    PP=Polyhedron(ieqs=ineqs,eqns=eqs)

    #Return the integral points
    return [PP.dim(),PP.integral_points()]

def Fill_Table_of_Lambdas(delta,ListP,tau : Tau,V : Representation,table): # delta is a sage vector in ZZ -  output : void
    """
    Fill the table with the possible entries for Lambda.
    Used for Kronecker
    """
    BL=[]
    for k in range(len(V.G)):
        for j in range(len(delta)):
            table[j,k]=Partitions(delta[j],max_length=tau.components[k][ListP[j][k]]).list()
    

def Fill_Table_of_Mus(delta,ListP, mtau : list[int],V : Representation,table): # delta is a sage vector in ZZ -  output : void
    """
    Fill the table with the possible entries for Lambda.
    Used for Fermion and Boson
    """
    #print('table',table)
    for k in range(len(mtau)):
        for j in range(len(delta)):
            #print('j,k:',j,k)
            table[j,k]=Partitions(ListP[j][k]*delta[j],max_length=mtau[k]).list()
    




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
                       
def Product_of_Tables(table_of_lists) : #TODO : déplacer dans utils
    """
    The entry is a table of lists. Produce the list of Tables obtained by picking (a product) one element in each list.
    """
    List_of_Tables=[]
    a,b=table_of_lists.shape
    
    for p in itertools.product(*table_of_lists.ravel()):
        T= np.empty(table_of_lists.shape, dtype=object)
        for i in range(a):
            for j in range(b):
                T[i,j]=p[i*b+j]
        List_of_Tables.append(T)
    return List_of_Tables

def Multiplicity_SV_tau(tau : Tau,chi : vector, V : Representation) -> int:

    #print('tau',tau)
    #print('chi',chi)
    ListP=tau.summands_Vtau(V)
    #print('ListP',ListP)
    Nu=chi2Numat(chi,tau.reduced.mult)  # Nu is a dominant weight for G^\tau as a table of partitions
    #print('Nu',Nu)
    w_Nu=fct_weights_of_Nu(Nu)
    #print('w_Nu',w_Nu)
    mult=0

    for delta in Enumerate_delta(ListP,w_Nu,V)[1]: # Run over the dela satisfying Condition 2
        #print('delta',delta)
        if V.type == 'kron' :
            table_Lambda = np.empty((len(ListP), len(V.G)), dtype=object)
            Fill_Table_of_Lambdas(delta,ListP,Tau(tau.reduced.mult),V,table_Lambda) # The list of tabulars of partitions satisfying Condition 1 and 3, for the given delta_j

            # Create the list of Lambda's from table_Lambda
            #List_of_Lambdas = [np.array(p, dtype=object).reshape(table_Lambda.shape) for p in itertools.product(*table_Lambda.ravel())]
            # reshape ne marche pas si p est une liste de listes à 1 éléments. Donc à la main:
        
            List_of_Lambdas=Product_of_Tables(table_Lambda)
            #a,b=table_Lambda.shape
    
            #for p in itertools.product(*table_Lambda.ravel()):
            #    T= np.empty(table_Lambda.shape, dtype=object)
            #    for i in range(a):
            #        for j in range(b):
            #            T[i,j]=p[i*b+j]
            #    List_of_Lambdas.append(T)
            
            for Lambda in List_of_Lambdas : 
                LR_coeff=LR_mat(Nu,Lambda,ListP)
                if LR_coeff!=0:
                    K_coeff=Kron_mat(Lambda)
                    mult+=LR_coeff*K_coeff
        else :
            s=len(tau.reduced.mult[0]) #[0] for the fist (and unique) bloc
            table_Mu = np.empty((len(ListP),s), dtype=object)
            #print('listP',len(ListP),ListP)
            Fill_Table_of_Mus(delta,ListP,tau.reduced.mult[0],V,table_Mu) # The list of tabulars of partitions satisfying Condition 1 and 3, for the given delta
            List_of_Mus=Product_of_Tables(table_Mu)
            for Mu in List_of_Mus:
                #print('delta',delta)
                #print('length ListP',ListP,len(ListP))
                #print('Mu',Mu.shape,Mu)
                #print('Nu',Nu)
                #print('Nu1O',Nu[1,0])
                lr=LR_mat_fb(Mu,Nu)
                if lr != 0 :
                  #if list(Mu[1,1])==[4,3,1] and list(Mu[1,0])==[5,3] and list(Mu[1,2])==[5,3] and list(Mu[0,1])==[3,3,3]: #TODO : A supprimer
                    #print('Mu ok:',Mu, lr)
                    #print('Mu :',Mu)    
                    Atot=1
                    sumAg=1
                    i=0
                    while sumAg !=0 and i<len(ListP):
                        sumAg=0
                        I=ListP[i]
                        if V.type == 'fermion':
                            Theta=[OurPartition(k*[1]) for k in I]
                        else :
                            Theta=[OurPartition([k]) for k in I]

                        #print('ici')    
                        for Lambda in itertools.product(OurPartition.all_for_integer(delta[i]),repeat=s) :
                            g=Kron_multi(Lambda)
                          #if list(Lambda[0])== [5,3] and list(Lambda[1])== [4,3,1] and list(Lambda[2])== [5,3] :#TODO : A supprimer
                          #  print('Lambda cou',Lambda,g)
                            if g!=0 :
                                a = Pleth_list(Lambda,Theta,Mu[i,:])
                                #if list(Lambda[1])== [4,3,1] and g == 2:
                                    #print('coucou')
                                    #print('a,g',a,g)
                                sumAg +=a*g 
                        i+=1
                    mult+=lr*sumAg    
                
        return(mult)            


    
def Is_Multiplicity_SV_tau_one(tau : Tau,chi : vector) -> bool: # TODO : version courte : calcule la multiplicité puis vérifie ==1. version longue et rapide : s'arrête dès que la multiplicité est >=2.

    d=tau.d
    delta=chi[0]
    ListP=tau.indices_in_tau_red_that_sum_to_zero
    Nu=chi2Numat(chi,tau.reduced.mult)  # Nu is a dominant weight for G^\tau

    w_Nu=fct_weights_of_Nu(Nu)
    mult=0
    Integral_Points=Enumerate_delta(delta,ListP,w_Nu,V)
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
