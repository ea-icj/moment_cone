import numpy as np
import itertools

from numpy.typing import NDArray

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
from .permutation import Permutation
from .kronecker import *

sym_f = SymmetricFunctions(QQ).s() 

def Search_Zero_a(Mu,ListP,Lambda,Vanishing_a): # Check if (Mu,Lambda) contains an already computer zero plethysm coefficient
    a,b=Mu.shape
    for i,j in itertools.product(range(a),range(b)):
        if (Mu[i,j],ListP[i][j],Lambda[i,j]) in Vanishing_a :
            return True
    return False    


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
    
# Kronecker coefficient of n-uplet of partitions using a multi-level cache
Kron_multi = KroneckerCoefficientMLCache()  


def all_partitions_of_max_length(n: int, l: Sequence[int], kro: KroneckerCoefficient) -> Iterable[tuple[tuple[Partition, ...], int]]:
    """
    All len(l)-uplets of partitions of n of non-zero Kronecker coefficient
    and so that len(p_i) <= l_i
    """
    # Case len(l)==1
    if len(l) == 1 and l[0] > 0:
        yield (Partition((n,)),),1
    elif len(l)==2:
        lm=min(l)
        for p in Partition.all_for_integer(n, lm):
            yield (p,p),1
    else :        
        # Sort by increasing maximal length (faster) and keep order
        permutation, sorted_l = zip(*sorted(enumerate(l), key=lambda l: l[1]))
        p_inverse = Permutation(permutation).inverse

        # All nuplets of partitions without the last length constraint
        head_product = itertools.product(*(Partition.all_for_integer(n, li) for li in sorted_l[:-1]))

        # Computing the product and yielding only the partitions of the decomposition
        # whose length respects the last constraint.
        for head in head_product:
            product = kro.product(head)
            for p, c in product.items():
                if len(p) <= sorted_l[-1]:
                    yield p_inverse(head + (p,)), c


def all_lambda_matrix(delta: Sequence[int], max_length: NDArray, kro: KroneckerCoefficient) -> Iterable[tuple[NDArray, int]]:
    """ All Lambda matrices form given weight vector and maximal length constraints
    
    Yield a matrix and the product of the Kronecker coefficient of each row.
    """
    N, s = max_length.shape
    assert N == len(delta)

    import numpy as np

    row_product = itertools.product(*(
        all_partitions_of_max_length(n, l, kro)
        for n, l in zip(delta, max_length)
    ))

    for rows_and_coeff in row_product:
        lambda_matrix = np.empty((N, s), dtype=object)
        lambda_matrix[:, :], coeffs = zip(*rows_and_coeff)
        yield lambda_matrix, prod(coeffs)
     

# TODO: dedicated file like for Kronecker?
class PlethysmCache:
    """ Cache for the Plethysm product """
    _cache: dict[tuple[Partition, Partition], dict[Partition, int]]
    _hit: int
    _miss: int

    def __init__(self):
        self._cache = dict()
        self._hit = 0
        self._miss = 0

    def __sym_f(self, p: Partition) -> Any:
        return sym_f(tuple(p))
    
    def product(self, lhs: Partition, rhs: Partition) -> dict[Partition, int]:
        try:
            result = self._cache[(lhs, rhs)]
            self._hit += 1
            return result
        except KeyError:
            pass

        self._miss += 1
        # Product using Sage
        product = self.__sym_f(lhs).plethysm(self.__sym_f(rhs))

        # Splitting decomposition
        result: dict[Partition, int] = dict()
        for monomial, coeff in product.monomial_coefficients().items():
            # Forcing int coefficients to remove Sage custom Integer type
            result[Partition(map(int, monomial))] = int(coeff)
        self._cache[(lhs, rhs)] = result
        return result
    
    def __call__(self, a: Partition, b: Partition, c: Partition) -> int:
        return self.product(a, b).get(c, 0)

    def __repr__(self) -> str:
        return f"PlethysmCache(#cache={len(self._cache)}, #hit={self._hit}, #miss={self._miss})"


plethysm_cache = PlethysmCache()

            
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
            Nu[i,k]=Partition(list(p))
            shift += mult
            
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
    #print('chi,tauredmult',chi,tau.reduced.mult)
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
                        to_be_deleted=[]
                        for i2,Lambda_tilde2 in  enumerate(reversed(List_of_Lambdas_plugged[i+1:])):
                            for j2 in range(len(ListP)) :
                                L2=Lambda_tilde2[0][j2,:]
                                if all(x ==y for x,y in zip(L2,L)) :
                                    to_be_deleted.append(i2)
                                    break
                        leng=len(List_of_Lambdas_plugged)
                        # TODO : ici on enregistre les effacements à faire puis on les fait. On peut sans doute effacer à mesure mais j'avais des erreurs. Sans doute pour des mauvaises raisons. 
                        for idx in to_be_deleted:
                            del List_of_Lambdas_plugged[leng-idx-1]
                        K_coeff=0
                        break # Unuseful to consider the other rows of Lambda_tilde
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
            
            List_of_Mus=Product_of_Tables(table_Mu)
            List_of_Mus_plugged=[]
            for Mu_tilde in List_of_Mus :
                LR=1
                Mu=np.empty((len(ListP),s), dtype=object)
                for j in range(s):
                    l=Mu_tilde[j,0]
                    LR*=l.mult
                    for i,la in enumerate(l.parts):
                            Mu[i,j]=la
                        
                List_of_Mus_plugged.append([Mu,LR])
            # Create the list of possible Lambda
            max_length=np.empty((len(ListP), s), dtype=int)
            for i in range(len(ListP)):
                for j,nb in enumerate(ListP[i]):
                   max_length[i,j]=Representation(LinGroup([tau.reduced.mult[0][j]]),V.type,nb_part=nb).dim
            #List_of_Lambdas,K=all_lambda_matrix(delta, max_length,1)
            
            # Runnig over the pairs Mu, Lambda #TODO : appelÃ©es mutilde et lambdatilde jeudi matin    
            Vanishing_a=set() # To remember the computed a that are zeros
            
            for [Mu,lr] in List_of_Mus_plugged: # lr is the multiplicity assocated to Mu
                #print('Next Mu',delta, max_length)
                
                for [Lambda,K] in  all_lambda_matrix(delta, max_length,Kron_multi): # K is the multiplicity assocated to Mu
                    #print('Next Lambda',Lambda)
                    if Search_Zero_a(Mu,ListP,Lambda,Vanishing_a): # Check if (Mu,Lambda) contains an already computer zero plethysm coefficient. In this case, we skip this pair.
                        break
                    A=1
                    for i,j in itertools.product(range(len(ListP)),range(s)):                                
                        if V.type == 'fermion':
                            theta=Partition(ListP[i][j]*[1])
                        else :
                            theta=Partition([ListP[i][j]])
                        #print('data plethysm',list(Lambda[i,j]),list(theta))
                        if len(Lambda[i,j])==max_length[i,j] :
                            #print('lambda can be reduced',Lambda[i,j],max_length[i,j])
                            la=Lambda[i,j].lambda_red(max_length[i,j])
                            if len(Lambda[i,j])==0:
                                shift=0
                            else :    
                                shift=theta[0]*Lambda[i,j][-1]
                            #print('La,theta,Mu',Lambda[i,j],theta,Mu[i,j])
                            if len(Mu[i,j])==0 or Mu[i,j][-1]<shift :
                                a=0
                            else :
                                mu=Partition([x-shift for x in Mu[i,j]])
                                a = plethysm_cache(la, theta, mu)
                                                            
                        else :
                            a = plethysm_cache(Lambda[i, j], theta, Mu[i, j])
                            #print('generic case',Lambda[i,j],theta,tau.reduced.mult[0][j])
                        
                        #print('a',a)
                        if a != 0 :
                            A*=a
                        else :
                            A=0
                            Vanishing_a.add((Mu[i,j],ListP[i][j],Lambda[i,j]))
                            break
                    mult+=lr*A*K             
                    if checkGreatEq2 and mult>1:
                        return(False)
                
    if checkGreatEq2:
        return(True)
    else :
        return(mult)            
