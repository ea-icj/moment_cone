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
from .kronecker import KroneckerCoefficientMLCache

sym_f = SymmetricFunctions(QQ).s()

class ListPartPlus:
    """
    A list of partitions with two more properties 
    - indices (a list of integers) and 
    - mult (an integer which is a multiplicity in representation theory.
    """

    def __init__(self, L : list[Partition], c : int,indices : Optional[list[int]]=None):
        """
        Initializes an instance of ListPartPlus.
        """
        self.parts=L
        self.mult=c
        self.indices=indices
        

    def __repr__(self) -> str:
        if self.indices != None :
            return 'Partitions: '+str(self.parts)+', Indices: '+str(self.indices)+', Multiplicity: '+str(self.mult)
        else :
            return 'Partitions: '+str(self.parts)+', Multiplicity: '+ str(self.mult)

    def __eq__(self,other : "ListPartPlus") -> bool:
        if all(l==m for l,m in zip(self.parts,other.parts)) and self.mult==other.mult and self.indices==other.indices :
            return True
        return False


def ListNonZeroLR(nu : Partition,delta : list[int],l:int):  
    s=len(delta)
    if not isinstance(nu, Partition):
        nu=Partition(nu)
    if s==1 :
        return([ListPartPlus([nu],1)])
    if s==2 :
        res=[]
        #print('nu',nu)
        nuc=nu.lambda_check(l).lambda_red(l)
        shift_max = delta[0]//l
        for a in range(shift_max+1):
            for lar in Partition.all_for_integer(delta[0]-a*l, max_length=l-1):
                #print('lambda reduit',lar)
                prod=lr.mult(nuc,lar)
                larn=Partition(lar)
                lamb1=Partition([larn[i]+a for i in range(l)])
                for mu in prod.keys():
                    if mu.length() <= l :
                        muc=Partition(mu).lambda_check(l)
                        b=(delta[1]-sum(muc))//l
                        #print(lar,lamb1,mu,muc,b)
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
                print(lr)
                if lr !=0:
                     res.append(ListPartPlus([Partition(l0),Partition(l1),Partition(l2)],lr))
    return(res)                 
        
def ListNonZeroLR_vtest2(nu : Partition,delta : list[int],l:int): # TODO : supprimer quand le reste marchera
    l=lr.coprod(nu)
    lt=[(x,l[x]) for x in l.keys() if len(x[0])<=l and len(x[1])<=l and sum(x[0])==delta[0]]
    return(lt)
    

# Kronecker coefficient of multiple partitions using a cache
Kron_multi: Callable[[Sequence[Partition]], int] = KroneckerCoefficientMLCache() 
        

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
    if shift!=len(chi):
       print(shift,len(chi), chi,mult_tau)
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
            table[j,k]=Partition.all_for_integer(delta[j],max_length=tau.components[k][ListP[j][k]]).list()
    

def Fill_Table_of_Mus(delta,ListP, mtau : list[int],V : Representation,table): # delta is a sage vector in ZZ -  output : void
    """
    Fill the table with the possible entries for Lambda.
    Used for Fermion and Boson
    """
    #print('table',table)
    for k in range(len(mtau)):
        for j in range(len(delta)):
            #print('j,k:',j,k)
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

def Multiplicity_SV_tau(tau : Tau,chi : vector, V : Representation, checkGreatEq2:bool=False) -> int:

    ListP=tau.summands_Vtau(V)
    Nu=chi2Numat(chi,tau.reduced.mult)  # Nu is a dominant weight for G^\tau as a table of partitions
    w_Nu=fct_weights_of_Nu(Nu)
    mult=0
    
    for delta in Enumerate_delta(ListP,w_Nu,V)[1]: # Run over the dela satisfying Condition 2
        print('delta',delta)
        if V.type == 'kron' :
            # Run over entries of Nu
            p,s = Nu.shape
            table_Lambda = np.empty((p,s), dtype=object)
            for k in range(s): # Two loops to run over the entries of Nu
                i=0
                while i < p and Nu[i,k] is not None:
                    Indices=[j for j,I in enumerate(listP) if ListP[j][k]==i]
                    max_length = Tau(tau.reduced.mult).components[k][i]
                    L=ListNonZeroLR(Nu[i,k],[delta[j] for j in Indices],max_length) # A list of ListPartPlus with Indicies=None
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
                for L in Lambda_tilde:
                    for k in range(s):
                        i=0
                        while i < p and L[i,k] is not None:
                            l=L[i,k]
                            LR*=l.mult
                            for j,la in zip(l.indices,l.parts):
                                Lambda[j,k]=la
                            i+=1
                List_of_Lambdas_plugged.append([Lambda,LR])
            i=0
            while i<len(List_of_Lambdas_plugge):
                Lambda_tilde = List_of_Lambdas_plugged[i]    
                K_coeff=Kron_mat(Lambda[0])
                if K !=0 :
                    mult+=Lambda[1]*K_coeff
                    if checkGreatEq2 and mult>1:
                        return(False)
                else:
                    ## Cancel the Lambda's having this row and hence giving K=0
                    #TODO
                    i=i
                i+=1    
        else :
            s=len(tau.reduced.mult[0]) #[0] for the first (and unique) bloc
            table_Mu = np.empty((s), dtype=object) # table_Mu[i] will be the list of possible columns i for Mu
            for i in range(s):
                max_length = Tau(tau.reduced.mult).components[0][i]
                table_Mu[i]=ListNonZeroLR(Nu[i,0],[ListP[j][i]*delta[j] for j in Indices],max_length)#coucou
            #Fill_Table_of_Mus(delta,ListP,tau.reduced.mult[0],V,table_Mu) # The list of tabulars of partitions satisfying Condition 1 and 3, for the given delta
            List_of_Mus=Product_of_Tables(table_Mu)
            for Mu in List_of_Mus:
                #print('delta',delta)
                #print('length ListP',ListP,len(ListP))
                #print('Mu',Mu.shape,Mu)
                #print('Nu',Nu)
                #print('Nu1O',Nu[1,0])
                lr=prod([l.mult for l in Mu]) # Each column of Mu has a multiplity we compute their product
                #if lr != 0 :
                  #if list(Mu[1,1])==[4,3,1] and list(Mu[1,0])==[5,3] and list(Mu[1,2])==[5,3] and list(Mu[0,1])==[3,3,3]: #TODO : A supprimer
                    #print('Mu ok:',Mu, lr)
                    #print('Mu :',Mu)    
                Atot=1
                sumAg=1
                i=0
                while sumAg !=0 and k<len(ListP):
                    sumAg=0
                    I=ListP[k]
                    if V.type == 'fermion':
                        Theta=[Partition(t*[1]) for t in I]
                    else :
                        Theta=[Partition([t]) for t in I]

                        #print('ici')    
                    for Lambda in itertools.product(Partition.all_for_integer(delta[i]),repeat=s) :
                        g=Kron_multi(Lambda)
                          #if list(Lambda[0])== [5,3] and list(Lambda[1])== [4,3,1] and list(Lambda[2])== [5,3] :#TODO : A supprimer
                          #  print('Lambda cou',Lambda,g)
                        if g!=0 :
                            a = Pleth_list(Lambda,Theta,Mu[i,:])
                                #if list(Lambda[1])== [4,3,1] and g == 2:
                                    #print('coucou')
                                    #print('a,g',a,g)
                            sumAg +=a*g
                            if checkGreatEq2 and mult+lr*sumAg >1 :
                                return(False)
                    i+=1
                mult+=lr*sumAg
                if checkGreatEq2 and mult>1:
                    return(False)
                
    if checkGreatEq2:
        return(True)
    else :
        return(mult)            
