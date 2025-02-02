from collections import Counter


def flatten(i,j,k,c=0):
    n=i+j+k
    T=(n+1)*(n+2)/2
    S = i*(2*n+3-i)/2+j
    return int(c*(T)+S)

#def unflatten(id:int,n:int,c : int=1):#TODO

def enumerate_vertices(n:int):
    for i in range(n + 1):
        for j in range(n - i + 1):
            k = n - i - j
            yield (i, j, k)

def rombus_ineq(n:int):
    T=int((n+1)*(n+2)/2)
    res=[]
    for (i,j,k) in enumerate_vertices(n-2) :
        L=T*[0]
        L[flatten(i,j+2,k)]=-1
        L[flatten(i+1,j,k+1)]=-1
        L[flatten(i+1,j+1,k)]=1
        L[flatten(i,j+1,k+1)]=1
        res.append([0]+L)
        L=T*[0]
        L[flatten(i+2,j,k)]=-1
        L[flatten(i,j+1,k+1)]=-1
        L[flatten(i+1,j+1,k)]=1
        L[flatten(i+1,j,k+1)]=1
        res.append([0]+L)
        L=T*[0]
        L[flatten(i,j,k+2)]=-1
        L[flatten(i+1,j+1,k)]=-1
        L[flatten(i+1,j,k+1)]=1
        L[flatten(i,j+1,k+1)]=1
        res.append([0]+L)
    return(res)    

def concave_ineq(n:int) -> list[list[int]]:
    """
    Inequalities asserting that the boundaries lambda and mu and concave and nondecreasing. That is correspond to partitions.
    """
    T=int((n+1)*(n+2)/2)
    for j in range(n):
        L=T*[0]
        L[flatten(0,j,n-j)]=-1
        L[flatten(0,j+1,n-j-1)]=1
        res=[[0]+L]
    for i in range(n):
        L=T*[0]
        L[flatten(i,n-i,0)]=-1
        L[flatten(i+1,n-i-1,0)]=1
        res.append([0]+L)        
    for j in range(1,n):
        L=T*[0]
        L[flatten(0,j,n-j)]=2
        L[flatten(0,j-1,n-j+1)]=-1
        L[flatten(0,j+1,n-j-1)]=-1
        res.append([0]+L)
    for i in range(1,n):
        L=T*[0]
        L[flatten(i,n-i,0)]=2
        L[flatten(i-1,n-i+1,0)]=-1
        L[flatten(i+1,n-i-1,0)]=-1
        res.append([0]+L)    
    return(res)
    
def glued_eq(n:int,c:int=0):
    """
    Gluing c+1 hives.
    """
    T=int((n+1)*(n+2)/2)
    list_ineq=[]
    list_eq=[[0,1]+(T*(c+1)-1)*[0]] # the (0,0,n) at 0 in the first hive
    # Rombus inequalities
    for l in range(c+1):
        for L in rombus_ineq(n):
            list_ineq.append([0]*(l*T+1)+L[1:]+[0]*T*(c-l))
    # gluing equalities
    for l in range(c):
        for i in range(n+1):
            L=[0]*(T*(c+1))
            # Equation i0(n-i) in l = 0i(n-i) in l+1
            L[flatten(i,0,n-i,l)]=1
            L[flatten(0,i,n-i,l+1)]=-1
            list_eq.append([0]+L)
    return [list_ineq,list_eq]      
    
    
def read_parts(H,n:int):
    """
    Read lambda,mu,nu from the hive H
    """
    # i=0
    Li=[H[flatten(0,j,n-j)] for j in range(n+1)]
    la=[b-a for a,b in itertools.pairwise(Li)]
    # j=0
    Lj=[H[flatten(i,0,n-i)] for i in range(n+1)]
    nu=[b-a  for a,b in itertools.pairwise(Lj)]
    # k=0
    Lk=[H[flatten(i,n-i,0)] for i in range(n+1)]
    mu=[b-a  for a,b in itertools.pairwise(Lk)]
    return([la,mu,nu])

def read_lambdai_from_glued (H,n:int,c:int):
    """
    Read lambda,mu,nu from the hive H
    """
    if c==0 : # Only one hyve
        # i=0
        Li=[H[flatten(0,j,n-j)] for j in range(n+1)]
        la=[b-a for a,b in itertools.pairwise(Li)]
        # k=0
        Lk=[H[flatten(i,n-i,0)] for i in range(n+1)]
        mu=[b-a  for a,b in itertools.pairwise(Lk)]
        return([la,mu])
    T=int((n+1)*(n+2)/2)*c
    # k=0
    Lk=[H[flatten(i,n-i,0)] for i in range(n+1)]
    la=[b-a  for a,b in itertools.pairwise(Lk)]
    return(read_lambdai_from_glued(H,n,c-1)+[la])
    

def init_nu(nu,n:int,c:int=0):
    """
    List of equalities putting nu on the boundary of the cth and last hive
    """
    nu+=(n-len(nu))*[0] #add zeros to get a partition of length n
    T=int((n+1)*(n+2)/2)
    res=[]
    s=0
    for i,t in enumerate(nu):
       L=T*(c+1)*[0]
       L[flatten(i,0,n-i,c)]=1
       res.append([-s]+L)
       s+=t
    L=(T*(c+1))*[0]
    L[flatten(n,0,0,c)]=1
    res.append([-s]+L)
    return(res)

def coprod_hives(nu,n:int ,c:int=0):
    T=int((n+1)*(n+2)/2)
    list_ineq,list_eq=glued_eq(n,c)
    print([len(L) for L in list_ineq])
    print([L for L in list_eq])
    list_concave_increasing=concave_ineq(n) 
    for i in range(c+1):
        for L in list_concave_increasing:
            L_shifted=(T*i+1)*[0]+L[1:]+T*(c-i)*[0] #The 1's to consider the constant term in the inequalities
            list_ineq.append(L_shifted)
    list_eq+=init_nu(nu,n,c)
    print([len(L) for L in list_ineq])
    print([len(L) for L in list_eq])
    #print(list_eq)
    #print('Ineq',list_ineq)
    PP=Polyhedron(ieqs=list_ineq,eqns=list_eq)
    res=[]
    #return(PP.integral_points_count())
    for H in PP.integral_points(): #TODO : ne marche que pour c=0 faire en général
        print('H',H)
        #la,mu,nu=read_parts(H,n) 
        lambdas=read_lambdai_from_glued (H,n,c)
        print('lambdas',lambdas)
        res.append(tuple([tuple(l) for l in lambdas])) # transformed in tuples of tuples to be hashable
        #res.append((tuple(la),tuple(mu)))
    print('len',len(res))    
    return(Counter(res))    

def init_lamu(la,mu,n:int): #TODO
    """
    List of equalities putting la and mu on the boundary 
    """
    mu+=n-len(mu)*[0] #add zeros to get  partitions of length n
    la+=n-len(la)*[0]
    T=(n+1)*(n+2)/2
    res=[]
    s=0
    for j,t in enumerate(la):
       L=T*[0]
       L[flatten(0,j,n-j)]=1
       res.append([-s]+L)
       s+=t
    L=T*[0]
    L[flatten(0,n,0)]=1
    res.append([-s]+L)
    for i,t in enumerate(mu):
       L=T*[0]
       L[flatten(i,n-i,0)]=1
       res.append([-s]+L)
       s+=t
    L=T*[0]
    L[flatten(n,0,0)]=1
    res.append([-s]+L)
    return(res)
    
class Hives_Ineq:

    def __init__(self, n: int,shpe : int=1):
        self.n=n

        

    
