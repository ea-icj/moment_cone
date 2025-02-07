import numpy as np
import itertools

from sage.all import SymmetricFunctions,ZZ,QQ,vector,matrix,Polyhedron # type: ignore
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


def ListNonZeroLR(nu : Partition,delta : list[int],l:int):
    s=len(delta)
    if not isinstance(nu, Partition):
        nu=Partition(nu)
    if s==0 :
        if len(nu) == 0 or nu[0]==0 :
            return([ListPartPlus([Partition([0])],1)])
        else :
            return([])
    if s==1 :
        return([ListPartPlus([nu],1)])
    if s==2 :
        res=[]
        nuc=nu.lambda_check(l).lambda_red(l)
        shift_max = delta[0]//l
        for a in range(shift_max+1):
            for lar in Partition.all_for_integer(delta[0]-a*l, max_length=l-1):
                prod=lr.mult(nuc,lar)
                larn=Partition(lar)
                lamb1=Partition([larn[i]+a for i in range(l)])
                for mu in prod.keys():
                    if mu.length() <= l :
                        muc=Partition(mu).lambda_check(l)
                        b=(delta[1]-sum(muc))//l
                        
                        if b+muc[l] >=0 :
                            lamb2=Partition([muc[i]+b for i in range(l)])
                            res.append(ListPartPlus([lamb1,lamb2],prod[mu]))
                    
        return(res)

    res=[]
    L2=ListNonZeroLR(nu,[delta[0],sum(nu)-delta[0]],l)
    # Extraction and uniquify of the mus TODO : améliorer la façon d'uniquifier
    List_mus1=[list_part_plus.parts[1] for list_part_plus in L2]
    List_mus=[]
    for mu in List_mus1 :
        if mu not in List_mus:
            List_mus.append(mu)
            
    for mu in List_mus :
        Ls_minus_one=ListNonZeroLR(mu,delta[1:],l)
        for list_lambda_queu in ListNonZeroLR(mu,delta[1:],l):
            for list_lambda_start in L2 :
                if list_lambda_start.parts[1]==mu:
                    res.append(ListPartPlus([list_lambda_start.parts[0]]+list_lambda_queu.parts,list_lambda_start.mult*list_lambda_queu.mult))
   
    ## Fusion of res by adding multiplicities TODO : ce que j'ai fait est moche et pas efficace pour zipper
    zipped_list=[]
    for obj in res:
        new_obj=ListPartPlus(obj.parts,0)
        if new_obj not in zipped_list :
            zipped_list.append(new_obj)
            
    for obj in res:
        for new_obj in zipped_list:
            if new_obj.parts == obj.parts:
                new_obj.mult+=obj.mult
                break       

    # Création de la nouvelle liste d'objets
    return(zipped_list)
    
def ListNonZeroLR_vtest(nu : Partition,delta : list[int],l:int): # TODO : supprimer quand le reste marchera
    s=len(delta)
    res=[]
    for l0 in Partition.all_for_integer(delta[0], max_length=l):
        for l1 in Partition.all_for_integer(delta[1], max_length=l):
            for l2 in Partition.all_for_integer(delta[2], max_length=l):
                lr=LR_multi([l0,l1,l2],nu)
                if lr !=0:
                     res.append(ListPartPlus([Partition(l0),Partition(l1),Partition(l2)],lr))
    return(res)                 
        
def ListNonZeroLR_vtest2(nu : Partition,delta : list[int],l:int): # TODO : supprimer quand le reste marchera
    l=lr.coprod(nu)
    lt=[(x,l[x]) for x in l.keys() if len(x[0])<=l and len(x[1])<=l and sum(x[0])==delta[0]]
    return(lt)
    

def Kron_multi(L) -> int :
    """
    L is a np.array or a list of Partitions of the same weight, length at least 2
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
            if Partition(list(monomial))==L[2]:
                return(coeff)
        return(0)    
    #### end unuseful
    
    # recursive part
    tot=0
    for monomial, coeff in product.monomial_coefficients().items():
        tot+=coeff*Kron_multi(list(L[2:])+[Partition(monomial)])
    return(tot)    
        

def LR_multi(L,nu) -> int: 
    
    """
    L is a list of partitions, nu is a partition
    return the multiple LR-coeffient
    """
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
        pl=sym_f(list(l)).plethysm(sym_f(list(t))).coefficient(list(n)) #TODO : utiliser un caché pour sym_f(list(l)).plethysm(sym_f(list(t)))
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
        for j in range(i,p):
             M[j,k]=-1        
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
            Nu[i,k]=Partition(p)
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


def Enumerate_delta(ListP : list[list[int]],sizenu,V : Representation,delta : int=0) : # output : of list of integer vectors. sizenu a partial matrix of int as list of columns. ListP ??
    """ 
    return the list of possible weight vectors delta. That is satisfying condition 2 of Proposition prop:multCVtau.
    delta is a vector of nonnegative integers of size len(ListP).
    """
    # Polyhedron in Z^l
    l=len(ListP)
    # Polyhedron in N^l
    ineqs=ineqs_in_Nl(l)
    eqs=[]
    if V.type == 'kron' : 
        # equalities corresponding to conditions 2
        for k in range(sizenu.ncols()): # Run over columns of Nu
            for i in range(sizenu.nrows()): # Run over the rows of the column
                if sizenu[i,k]>=0:
                    v=[-sizenu[i,k]]
                    for j in range(l):
                        if ListP[j][k]==i:
                            v.append(1)
                        else:
                            v.append(0)
                    eqs.append(tuple(v))
    else:
        
        for j in range(sizenu.nrows()):
            eqs.append(tuple([-sizenu[j,0]]+[I[j] for I in ListP]))            
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
            table[j,k]=Partition.all_for_integer(delta[j],max_length=tau.components[k][ListP[j][k]]).list()
    

def Fill_Table_of_Mus(delta,ListP, mtau : list[int],V : Representation,table): # delta is a sage vector in ZZ -  output : void
    """
    Fill the table with the possible entries for Lambda.
    Used for Fermion and Boson
    """
    
    for k in range(len(mtau)):
        for j in range(len(delta)):
            table[j,k]=Partition.all_for_integer(ListP[j][k]*delta[j],max_length=mtau[k]).list()
    




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
    The entry is a bi-dimesional table of ListPartPlus. Produce the list of Tables obtained by picking (a product) one element in each list.
    Let entries None in place. 
    """
    List_of_Tables=[]
    a,b=table_of_lists.shape

    Inone=[]
    for i in range(a):
        for j in range(b):
            if table_of_lists[i,j]==None:
                table_of_lists[i,j]=[0]
                Inone.append([i,j])
                
    for p in itertools.product(*table_of_lists.ravel()):
        T= np.empty(table_of_lists.shape, dtype=object)
        for i in range(a):
            for j in range(b):
                T[i,j]=p[i*b+j]
        for [i,j] in Inone:
            T[i,j]=None
        List_of_Tables.append(T)
    return List_of_Tables

def Multiplicity_SV_tau(tau : Tau,chi : vector, V : Representation, checkGreatEq2:bool=False) -> int:
    from math import comb # equivalent to sage.binomial

    ListP=tau.summands_Vtau(V)
    #print('tau',tau)
    #print('ListP',ListP)
    #print('chi,tauredmult',chi,tau,tau.reduced.mult)
    Nu=chi2Numat(chi,tau.reduced.mult)  # Nu is a dominant weight for G^\tau as a table of partitions
    w_Nu=fct_weights_of_Nu(Nu)
    #print('Nu',Nu)
    #print('sizeNu',Nu.shape)
    #print('w_Nu',w_Nu)
    mult=0
    Delta=Enumerate_delta(ListP,w_Nu,V)
    if checkGreatEq2 and tau.is_dom_reg : # In this case we only need to check the dimension of the polyhedron of delta's
        return Delta[0]==0
    
    for delta in Delta[1]: # Run over the dela satisfying Condition 2
        #print('delta',delta)
        if V.type == 'kron' :
            # Run over entries of Nu
            p,s = Nu.shape
            table_Lambda = np.empty((p,s), dtype=object)
            for k in range(s): # Two loops to run over the entries of Nu
                i=0
                while i < p and Nu[i,k] is not None:
                    Indices=[j for j,I in enumerate(ListP) if ListP[j][k]==i]
                    max_length = tau.reduced.mult[k][i]
                    L=ListNonZeroLR(Nu[i,k],[delta[j] for j in Indices],max_length) # A list of ListPartPlus with Indices=None
                    ## Add indices
                    for l in L :
                        l.indices=Indices
                    table_Lambda[i,k] = L
                    i+=1

            List_of_Lambdas=Product_of_Tables(table_Lambda)
            List_of_Lambdas_plugged=[]
            for Lambda_tilde in List_of_Lambdas :
                LR=1
                Lambda=np.empty((len(ListP),s), dtype=object)
                
                for k in range(s):
                    i=0
                    while i < p and Lambda_tilde[i,k] is not None:
                        l=Lambda_tilde[i,k]
                        LR*=l.mult
                        for j,la in zip(l.indices,l.parts):
                            Lambda[j,k]=la
                        i+=1
                List_of_Lambdas_plugged.append([Lambda,LR])
                
            i=0
            while i<len(List_of_Lambdas_plugged):
                Lambda_tilde = List_of_Lambdas_plugged[i]
                K_coeff=1
                for j in range(len(ListP)) :
                    L=Lambda_tilde[0][j,:]
                    K=Kron_multi(L)
                    if K !=0 :
                        K_coeff*=K                                                
                    else:
                        ## Cancel the Lambda's having this row and hence giving K=0
                        #TODO : peut-on faire mieux que ça ?
                        #last=len(List_of_Lambdas_plugged)-1
                        to_be_deleted=[]
                        for i2,Lambda_tilde2 in  enumerate(reversed(List_of_Lambdas_plugged[i+1:])):
                            for j2 in range(len(ListP)) :
                                L2=Lambda_tilde2[0][j2,:]
                                if all(x ==y for x,y in zip(L2,L)) :
                                    #print('L2,L',L2,L)
                                    to_be_deleted.append(i2)
                                    #print('supprim')
                                    #print('i,i2,len',i,i2,len(List_of_Lambdas_plugged),len(List_of_Lambdas_plugged[i+1:]))
                                    #del List_of_Lambdas_plugged[-(i2+1)] #Remove Lambda_tilde2 from List_of_Lambdas_plugged
                                    #TODO : la ligne ci-dessus compile mais introduit une erreur. A corriger.
                                #else:
                                #    print('not supprim')
                                    break
                        #print([x[1] for x in List_of_Lambdas_plugged])
                        #print('to be deleted',to_be_deleted)
                        leng=len(List_of_Lambdas_plugged)        
                        for idx in to_be_deleted:
                            del List_of_Lambdas_plugged[leng-idx-1]
                        #print([x[1] for x in List_of_Lambdas_plugged])  
                        K_coeff=0
                        break
                mult+=Lambda_tilde[1]*K_coeff    
                if checkGreatEq2 and mult>1:
                    return(False)        
                i+=1    
        else :
            s=len(tau.reduced.mult[0]) #[0] for the first (and unique) bloc
            table_Mu = np.empty((s,1), dtype=object) # table_Mu[i] will be the list of possible columns i for Mu
            for j in range(s):
                max_length = tau.reduced.mult[0][j]
                table_Mu[j,0]=ListNonZeroLR(Nu[j,0],[ListP[i][j]*delta[i] for i in range(len(ListP))],max_length)
            #Fill_Table_of_Mus(delta,ListP,tau.reduced.mult[0],V,table_Mu) # The list of tabulars of partitions satisfying Condition 1 and 3, for the given delta
            List_of_Mus=Product_of_Tables(table_Mu)
            List_of_Mus_plugged=[]
            for Mu_tilde in List_of_Mus :
                #print('Lambda_tilde',Lambda_tilde)
                LR=1
                Mu=np.empty((len(ListP),s), dtype=object)
                for j in range(s):
                    l=Mu_tilde[j,0]
                    LR*=l.mult
                    for i,la in enumerate(l.parts):
                            Mu[i,j]=la
                        
                List_of_Mus_plugged.append([Mu,LR])
                
            for [Mu,lr] in List_of_Mus_plugged:
                Atot=1
                sumAg=1
                k=0
                while sumAg !=0 and k<len(ListP):
                    sumAg=0
                    I=ListP[k]
                    if V.type == 'fermion':
                        Theta=[Partition(t*[1]) for t in I]
                    else :
                        Theta=[Partition([t]) for t in I]

                    for Lambda in itertools.product(*[Partition.all_for_integer(delta[k],max_length=comb(tau.reduced.mult[0][j],I[j])) for j in range(s)]) :
                        g=Kron_multi(Lambda)
                        if g!=0 :
                            a = Pleth_list(Lambda,Theta,Mu[k,:])
                            sumAg +=a*g
                            if checkGreatEq2 and mult+lr*sumAg >1 :
                                return(False)
                    k+=1
                mult+=lr*sumAg
                if checkGreatEq2 and mult>1:
                    return(False)
                
    if checkGreatEq2:
        return(True)
    else :
        return(mult)            
