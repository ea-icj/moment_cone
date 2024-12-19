from cone.typing import *
from cone.dimension import *
from cone.vector_chooser import *
from cone.hyperplane_candidates import *
from cone.list_of_W import *
from cone.inequality import *
from cone.tau import *
from cone.ramification import *
from cone.Normaliz2python import *
from cone.dim_stabilizer_K import *

import sys

####
d0=Dimension([4,4,4])
tpi_method: Method ='symbolic'
#tpi_method='probabilistic'
ram_schub_method: Method = 'symbolic'
ram0_method: Method = 'probabilistic'
####



## Checking if the cone has the expected dimension

Ms=Lie_action_as_matrices_V(d0)
MsR=[mat_C_to_R(M) for M in Ms]
if dim_gen_stab_of_K(MsR)>len(d0)-1: # Check that the dim is computed in U_n(C)^s without the isolated S^1
    print('The moment cone has codimension bigger that the length of d. Namely',dim_gen_stab_of_K(MsR))
    print('The program does not work in this case')
    sys.exit()
else:
    print('The cone has the expected dimension.')

## Generate the list of candidates for tau

#Candidates_for_tau=find_hyperplanes_mod_sym_dim(d0,d0.dimU) # This is the function for regular ops (todo : include this info in the name) - To be changed.
print('Step 1, looking for a first list of dominant 1-PS whose kernel is supported at hyperplanes of weights.')

Candidates_for_tau_1ps =find_1PS_mod_sym_dim(d0)

print(len(Candidates_for_tau_1ps), ' dominant 1-PS selected in Step 1')
for tau in Candidates_for_tau_1ps:
    print(tau)
Candidates_for_tau=unique_modulo_symmetry_list_of_tau(Candidates_for_tau_1ps) # todo : inutile car déjà fait dans find_1PS ?
print(len(Candidates_for_tau), ' dominant 1-PS selected in Step 1 After Unicity')

for tau in Candidates_for_tau:
    print(tau)
    
# Filter 1: submodule condition

print('Step 2, Checking submodule condition')
Candidates_for_tau1=[tau for tau in Candidates_for_tau if tau.is_sub_module]
print(len(Candidates_for_tau1), ' dominant 1-PS selected in Step 2')

# Filter 2: stabilizer condition
print('Step 3, Stabilizer condition')


### Avec le nouveau dimStab
Candidates_for_tau2=[]
for tau in Candidates_for_tau1:
    if  tau.is_dom_reg :
        Candidates_for_tau2.append(tau)
    else: 
        Ms_tau=Lie_action_as_matrices_Vtau(tau,Ms)
        Ms_tauR=[mat_C_to_R(M) for M in Ms_tau]
        if dim_gen_stab_of_K(Ms_tauR)==len(d0):
            Candidates_for_tau2.append(tau)    
print(len(Candidates_for_tau2), ' dominant 1-PS selected in Step 3')


for tau in Candidates_for_tau2:
    print(tau)


## Generate the list of candidates for the inequalites (pairs tau,w)
## Here w has to belong to P^tau and U(w) is tau-isomorphic to V(tau>0)

print('Step 4, computing action of W, in order to get a first list of inequalities containing all the expected ones')
Candidates_for_Ineq=[]
for tau in Candidates_for_tau2 :
    Lw=ListWs_Mod(tau)
    Candidates_for_Ineq+=[Inequality(tau,w) for w in Lw] # Fait-on un dictionnaire tau : Lw ??
print(len(Candidates_for_Ineq), ' inequalities selected in Step 4')

for ineq in Candidates_for_Ineq :
    print(ineq)
    
# Filter 0: Unicity modulo sym(d)
print('Step 5, Reduction modulo symmetries of the dimension vector')
Candidates_for_Ineq1=unique_modulo_symmetry_list_of_ineq(Candidates_for_Ineq)
print(len(Candidates_for_Ineq1), ' inequalities selected in Step 5')

for ineq in Candidates_for_Ineq1 :
    print(ineq.tau)

# Filter 1: pi is dominant
print('Step 6, checking dominancy of the map pi')
Dominant_Ineq=[ineq for ineq in Candidates_for_Ineq1 if Check_Rank_Tpi(ineq,tpi_method)] 
print(len(Dominant_Ineq), ' inequalities selected in Step 6')

for ineq in Dominant_Ineq :
    print(ineq.tau)

# Filter 2: pi is birational (ramification divisor contracted)
print('Step 7, checking birationality (ramification divisor contracted) of the map pi')
Birational_Ineq=[ineq for ineq in Dominant_Ineq if Is_Ram_contracted(ineq,ram_schub_method,ram0_method)]
print(len(Birational_Ineq), ' inequalities selected in Step 7')

#path="/home/bm29130h/Documents/Recherche/Ressources_autres/GDT/Machine Learning/calculs Kron/2 oct/"
#reference=[Inequality.from_tau(tau) for tau in convert_file_Nout2pyth(path,d0)]
#dictionary_list_lengths(compare_ineq_candidates_reference_mod_sym_dim(Candidates_for_Ineq1,reference))
#test fails following fusion 17dec 12h11

#unique_reference=unique_modulo_symmetry_list_of_ineq(reference)
#dictionary_list_lengths(compare_tau_candidates_reference_mod_sym_dim(Candidates_for_tau,[ineq.tau for ineq in reference]))

for ineq in Birational_Ineq :
    print(ineq)
