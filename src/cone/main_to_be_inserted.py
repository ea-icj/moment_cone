from cone.typing import *
from cone.dimension import *
from cone.cone_dimension import *
from cone.hyperplane_candidates import *
from cone.listW_temp import *
from cone.inequality import *
from cone.tau import *
from cone.ramification import *


####
d0=Dimension([4,4,4])
#stabilizer_method='symbolic'
stabilizer_method='probabilistic'
tpi_method='symbolic'
ram_schub_method='probabilistic'
ram0_method='probabilistic'
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
print('Step 1, looking for a first list of dominant 1-PS whose kernel is supported at hyperplanes of weights.')
Candidates_for_tau=find_1PS_mod_sym_dim(d0)   
print(len(Candidates_for_tau), ' dominant 1-PS selected in Step 1')
# Filter 1: submodule condition

print('Step 2, Checking submodule condition')
Candidates_for_tau1=[tau for tau in Candidates_for_tau if tau.is_sub_module]
print(len(Candidates_for_tau1), ' dominant 1-PS selected in Step 2')

# Filter 2: stabilizer condition
print('Step 3, Stabilizer condition')
Candidates_for_tau2=[tau for tau in Candidates_for_tau1 if dim_of_stabilizer_in_K_tau(tau,stabilizer_method)==len(d0)]
print(len(Candidates_for_tau2), ' dominant 1-PS selected in Step 3')

## Generate the list of candidates for the inequalites (pairs tau,w)
## Here w has to belong to P^tau and U(w) is tau-isomorphic to V(tau>0)
print('Step 4, computing action of W, in order to get (possibly non-dominant) 1-PS, yielding a redundant list of inequalities for the cone')
Candidates_for_Ineq=[]
for tau in Candidates_for_tau2 :
    Lw=ListWs_Mod(tau)
    Candidates_for_Ineq+=[Inequality(tau,w) for w in Lw] # Fait-on un dictionnaire tau : Lw ??
print(len(Candidates_for_Ineq), ' inequalities selected in Step 4')

# Filter 0: Unicity modulo sym(d)
print('Step 5, Reduction modulo symmetries of the dimension vector')
Candidates_for_Ineq1=unique_modulo_symmetry_list_of_ineq(Candidates_for_Ineq)
print(len(Candidates_for_Ineq1), ' inequalities selected in Step 5')

# Filter 1: pi is dominant
print('Step 6, checking dominance of map pi')
Dominant_Ineq=[ineq for ineq in Candidates_for_Ineq1 if Check_Rank_Tpi(ineq,tpi_method)] 
print(len(Dominant_Ineq), ' inequalities selected in Step 6')

# Filter 2: pi is dominant
print('Step 7, checking dominance of map pi (again)')
Birational_Ineq=[ineq for ineq in Dominant_Ineq if Is_Ram_contracted(ineq,ram_schub_method,ram0_method)]
print(len(Birational_Ineq), ' inequalities selected in Step 7')


