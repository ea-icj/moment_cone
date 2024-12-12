from src.cone.typing import *
from src.cone.dimension import *
from src.cone.cone_dimension import *
from src.cone.hyperplane_candidates import *



####
d0=Dimension([3,3,3])
#stabilizer_method='symbolic'
stabilizer_method='probabilistic'
tpi_method='symbolic'
####




## Checking if the cone has the expected dimension

tau_triv=Tau.from_flatten((d0.sum+1)*[0],d0) #The trivial one parameter subgroup
if dim_of_stabilizer_in_K_tau(tau_triv,stabilizer_method)>len(d0)-1: # Check that the dim is computed in U_n(C)^s without the isolated S^1
    print('The moment cone has codimension bigger that the length of d. Namely', dim_of_stabilizer_in_K_tau(tau_triv,stabilizer_method))
    print('The program does not work in this case')
    sys.exit()
else:
    print('The cone has the expected dimension.')

## Generate the list of candidates for tau

#Candidates_for_tau=find_hyperplanes_mod_sym_dim(d0,d0.dimU) # This is the function for regular ops (todo : include this info in the name) - To be changed.
   # Cette est-elle déjà filrée modulo répétition et symétrie. Sinon faire un Filter 0.
   # En fait donne les listes de poids pas les tau. A modifier donc
Candidates_for_tau=find_1PS_mod_sym_dim(d0)   
# Filter 1: submodule condition



Candidates_for_tau1=[tau for tau in Candidates_for_tau if tau.Check_sub_module()]

# Filter 2: stabilizer condition

Candidates_for_tau2=[tau for tau in Candidates_for_tau1 if dim_of_stabilizer_in_K_tau(tau,stabilizer_method)==len(d0)]

## Generate the list of candidates for the inequalites (pairs tau,w)
## Here w has to belong to P^tau and U(w) is tau-isomorphic to V(tau>0)

Candidates_for_Ineq=[]
for tau in Candidates_for_tau2 :
    Lw=Candidates_for_tau2(tau)
    Candidates_for_Ineq+=[Inequality(tau,w) for w in Lw] # Fait-on un dictionnaire tau : Lw ??

# Filter 0: Unicity modulo sym(d)

# Filter 1: pi is dominant

Dominant_Ineq=[ineq for ineq in Candidates_for_Ineq if Check_Rank_Tpi(ineq,tpi_method)] 
