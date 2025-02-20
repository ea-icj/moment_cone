from .typing import *
from .linear_group import *
from .representation import *
from .stabK import *
from .tau import find_1PS
from .inequality import *
from .list_of_W import *
from .linear_triangular import *
from .bkr import *
from .ramification import *
from .Groebner import *
from .export import *
from .task import *


G = LinearGroup([3, 3, 3, 1])
V = KroneckerRepresentation(G)
#G = LinearGroup([8])
#V = FermionRepresentation(G, particle_cnt=3)
#G = LinearGroup([4])
#V = BosonRepresentation(G, particle_cnt=2)

print("We are computing the inequalities for the moment cone of",V)
#tpi_method: Method ='symbolic'
tpi_method: Method ='probabilistic'
#ram_schub_method: Method = 'symbolic'
ram_schub_method: Method = 'probabilistic'
ram0_method: Method = 'probabilistic'
#grobner_method: Method = 'symbolic'
#grobner_lim=3 #limit of time for each Grobner computation
grobner_method: Method = 'probabilistic'
grobner_lim=2

## Checking if the cone has the expected dimension

Ms=V.actionK

MsR = [mat_C_to_R(M) for M in Ms.values()]
        
#MsR=[mat_C_to_R(M) for M in Ms]
if dim_gen_stab_of_K(MsR)>G.rank-V.dim_cone: # Check that the dim is computed in U_n(C)^s without the isolated S^1
    print('The general stabilizer of K in V is too big. Namely of dimension',dim_gen_stab_of_K(MsR))
    print('The program does not work in this case')
    import sys
    sys.exit()
else:
    print('The cone has the expected dimension.')

## Generate the list of candidates for tau

print('Step 1, looking for a first list of dominant 1-PS whose kernel is supported at hyperplanes of weights.')

Candidates_for_tau = find_1PS(V)

print(len(Candidates_for_tau), ' dominant 1-PS with hyperplane spanned by weights or no too many positive weights')


# Filter 1: submodule condition

print('Step 2, Checking submodule condition')
Candidates_for_tau1=[tau for tau in Candidates_for_tau if tau.is_sub_module(V)]
print(len(Candidates_for_tau1), ' dominant 1-PS satisfying the submodule condition')


# Filter 2: stabilizer condition
print('Step 3, Stabilizer condition')

Candidates_for_tau2=[]
for tau in Candidates_for_tau1:
    if  tau.is_dom_reg :
        Candidates_for_tau2.append(tau)
    else: 
        Ms_tau=Lie_action_as_matrices_Vtau(tau,Ms,V)
        Ms_tauR=[mat_C_to_R(M) for M in Ms_tau.values()]
        
        if dim_gen_stab_of_K(Ms_tauR)==G.rank-V.dim_cone+1:
            Candidates_for_tau2.append(tau)    
print(len(Candidates_for_tau2), ' dominant 1-PS satisfying the stabilizer condition')

#print(Candidates_for_tau2)

## Generate the list of candidates for the inequalites (pairs tau,w)
## Here w has to belong to P^tau and U(w) is tau-isomorphic to V(tau>0)

print('Step 4, for each tau, computation the w with compatible tau-modules')
Candidates_for_Ineq=[]
for tau in Candidates_for_tau2 :
    Lw=ListWs_Mod(tau,V)
    Candidates_for_Ineq+=[Inequality(tau,w) for w in Lw] # TODO : Fait-on un dictionnaire tau -> liste de w ?
print(len(Candidates_for_Ineq), ' inequalities selected in Step 4')

print('Step 5, Reduction modulo symmetries of the dimension vector')
#for ineq in Candidates_for_Ineq:
#    print('ineq',ineq)
#    print('sorted',ineq.sort_mod_sym_dim)
    
Candidates_for_Ineq1=list(set(ineq.sort_mod_sym_dim for ineq in Candidates_for_Ineq))
#unique_modulo_symmetry_list_of_ineq()
print(len(Candidates_for_Ineq1), ' inequalities selected in Step 5')

# Filter 1: pi is dominant
print('Step 6, checking dominancy of the map pi')
Dominant_Ineq=[ineq for ineq in Candidates_for_Ineq1 if Check_Rank_Tpi(ineq,V,tpi_method)] 
print(len(Dominant_Ineq), ' inequalities selected in Step 6')

# Filter 2: Linear Triangular
print('Step 7, checking Linear Triangular')
Ineq_Triang=[]
Dominant_Ineq_filteredLT=[]
for ineq in Dominant_Ineq:
    if is_linear_triangular(V,ineq.tau,[alpha for alpha in ineq.inversions]) :
        Ineq_Triang.append(ineq)
    else:
        Dominant_Ineq_filteredLT.append(ineq)
print('There are',len(Ineq_Triang),'Linear Triangular inequalities')

############################ possible structure for a timeout for BKR. Not working yet
#print('alternative 8 Step, introducing timeouts for BKR')
#List_BKR=[]
#List_BKR_inc=[]
#for ineq in Dominant_Ineq_filteredLT:
#     try:
#        with timeout(1, no_raise=False):
#            chi=ineq.weight_det(V)
#            if list(ineq.inversions)==[] or Multiplicity_SV_tau(ineq.tau,chi,V,True):    
#                List_BKR.append(ineq)
#                print(ineq, "selected in List_BKR (True)")
#            else:
#                print(ineq,"not selected in List_BKR (False)")
#     except TimeOutException:
#         List_BKR_inc.append(ineq)
#         print(ineq, "inconclusive computation")
#
#print('There are',len(List_BKR), 'selected inequalities in List_BKR and ',len(List_BKR_inc), 'inconclusive inequalities')
 
############################
       
## Filter 3: BKR condition
if isinstance(V,ParticleRepresentation) and V.G[0]>=8:
    print('Step 8 skipped (condition BKR) for particle representation in rank >=8')
    List_BKR=Dominant_Ineq_filteredLT
else:
    print('Step 8, checking if BKR condition is fullfilled')
    List_BKR=[]
    for ineq in Dominant_Ineq_filteredLT :
        chi=ineq.weight_det(V)
        if list(ineq.inversions)==[] or Multiplicity_SV_tau(ineq.tau,chi,V,True):    
            List_BKR.append(ineq)
    print('The BKR filter has eliminated',len(Dominant_Ineq_filteredLT)-len(List_BKR),'inequalities')
    print("Final state of the Kronecker cache:", Kron_multi)
    print("Final state of the Plethysm cache: ", plethysm_cache)

        
#print('The BKR filter has eliminated',len(Dominant_Ineq_filteredLT)-len(List_BKR),'inequalities')
#print("Final state of the Kronecker cache:", Kron_multi)
#print("Final state of the Plethysm cache: ", plethysm_cache)

# Filter 4: pi is birational (ramification divisor contracted)
print('Step 9, checking birationality (ramification divisor contracted) of the map pi')
#Birational_Ineq=[ineq for ineq in Dominant_Ineq if Is_Ram_contracted(ineq,V,ram_schub_method,ram0_method)]
Birational_Ineq=[ineq for ineq in List_BKR if Is_Ram_contracted(ineq,V,ram_schub_method,ram0_method)]+Ineq_Triang
#Birational_Ineq=[ineq for ineq in Dominant_Ineq_filteredLT if Is_Ram_contracted(ineq,V,ram_schub_method,ram0_method)]
print(len(Birational_Ineq), ' inequalities selected in Step 9 in','seconds')


#TODO adapt PolynomialRingForWeights (especially variable_name) to roots in order to allow Grobner computations
############################ possible structure for a timeout for Groebner. Not working yet
print('alternative 8-9 Step, checking birationality via Grobner')
for ineq in List_BKR:
     try:
        with timeout(1, no_raise=False):
             print(is_fiber_singleton(V,ineq,"probabilistic"))
     except TimeOutException:
         print("False?")
########################### or its alternative
#print('alternative 8-9 Step, checking birationality via Grobner')
#Grobner_output=Grobner_List_Test(Dominant_Ineq_filteredLT,grobner_lim,V,grobner_method)
#True_Ineq=Grobner_output[0]+Ineq_Triang
#Dominant_Ineq_filteredGrobner=Grobner_output[1]
#print(len(True_Ineq), 'true inequalities after Grobner; presumably the only ones but', len(Dominant_Ineq_filteredGrobner), 'inequalities where Grobner was inconclusive')
###########################

         


#for ineq in Ineq_Triang :
#    if ineq not in Birational_Ineq:
#        print('Tri Lin pas birational:',ineq,list(ineq.inversions),ineq.tau.positive_weights(V))
#print('Triang')
#for ineq in Ineq_Triang :
#    print(ineq)

print('Birational_Ineq')

for ineq in Birational_Ineq :
    print(ineq)
#    chi=ineq.weight_det(V)
#    print('mult',Multiplicity_SV_tau(ineq.tau,chi,V))

#for ineq in List_BKR:
#     with timeout(3):
#         print(is_fiber_singleton(V,ineq,"probabilistic"))




#exports possibles:
export_normaliz(V,full_under_symmetry_list_of_ineq(Birational_Ineq),add_equations="all",add_dominance="all")
export_latex(V,Birational_Ineq,sgn=-1)
export_python(V,Birational_Ineq)


#comparaison possible:
#from src.All_Cases_Init.reference_datas.comparisons import *
#compare_to_reference(Birational_Ineq,V)




