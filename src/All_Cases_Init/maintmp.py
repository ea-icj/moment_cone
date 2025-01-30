
G = LinGroup([3,3,3,1])
V = Representation(G,'kron')
G = LinGroup([7])
V = Representation(G,'fermion',nb_part=3)

tpi_method: Method ='symbolic'
tpi_method='probabilistic'
ram_schub_method: Method = 'symbolic'
ram_schub_method: Method = 'probabilistic'
ram0_method: Method = 'probabilistic'

## Checking if the cone has the expected dimension

Ms=V.actionK

MsR = [mat_C_to_R(M) for M in Ms.values()]
        
#MsR=[mat_C_to_R(M) for M in Ms]
if dim_gen_stab_of_K(MsR)>G.rank-V.dim_cone: # Check that the dim is computed in U_n(C)^s without the isolated S^1
    print('The general stabilizer of K in V is too big. Namely of dimension',dim_gen_stab_of_K(MsR))
    print('The program does not work in this case')
    sys.exit()
else:
    print('The cone has the expected dimension.')

## Generate the list of candidates for tau

print('Step 1, looking for a first list of dominant 1-PS whose kernel is supported at hyperplanes of weights.')

Candidates_for_tau = find_1PS(V)

#print(Candidates_for_tau)
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

## Generate the list of candidates for the inequalites (pairs tau,w)
## Here w has to belong to P^tau and U(w) is tau-isomorphic to V(tau>0)

print('Step 4, for each tau, compution the w with compatible tau-modules')
Candidates_for_Ineq=[]
for tau in Candidates_for_tau2 :
    Lw=ListWs_Mod(tau,V)
    Candidates_for_Ineq+=[Inequality(tau,w) for w in Lw] # TODO : Fait-on un dictionnaire tau -> liste de w ?
print(len(Candidates_for_Ineq), ' inequalities selected in Step 4')

print('Step 5, Reduction modulo symmetries of the dimension vector')
Candidates_for_Ineq1=list(set(ineq.sort_mod_sym_dim for ineq in Candidates_for_Ineq))
#unique_modulo_symmetry_list_of_ineq()
print(len(Candidates_for_Ineq1), ' inequalities selected in Step 5')

# Filter 1: pi is dominant
print('Step 6, checking dominancy of the map pi')
Dominant_Ineq=[ineq for ineq in Candidates_for_Ineq1 if Check_Rank_Tpi(ineq,V,tpi_method)] 
print(len(Dominant_Ineq), ' inequalities selected in Step 6')

# Filter 3: BKR condition
print('Step 8, checking if BKR condition is fullfilled')
List_BKR=[]
for ineq in Dominant_Ineq :
    #print('ineq',ineq)
    #print('pos weights:',[chi for chi in ineq.tau.positive_weights(V)])
    #print('pos roots:',[chi for chi in ineq.inversions])
    chi=ineq.weight_det(V)
    print('Mult Vtau', Multiplicity_SV_tau(ineq.tau,chi,V)) 
    #if list(ineq.inversions)==[] or Is_Multiplicity_SV_tau_one(ineq.tau,chi):
    #    List_BKR.append(ineq)
        
# Filter 4: pi is birational (ramification divisor contracted)
print('Step 9, checking birationality (ramification divisor contracted) of the map pi')
Birational_Ineq=[ineq for ineq in Dominant_Ineq if Is_Ram_contracted(ineq,V,ram_schub_method,ram0_method)]
print(len(Birational_Ineq), ' inequalities selected in Step 9 in','seconds')








###### TODO : TESTS A SUPPRIMMER CI-APRES #########

for ineq in Birational_Ineq :
    chi=ineq.weight_det(V)
    print('Mult Vtau', Multiplicity_SV_tau(ineq.tau,chi,V))
    
tau=Tau.from_flatten([2,1,0,2,1,0,2,1,0,0],G)
for chi in V.all_weights:
    for chi2 in V.all_weights:
        if chi2.leq(chi) and chi != chi2:
            print(tau.dot_weight(chi)-tau.dot_weight(chi2))
#for ineq in Birational_Ineq :
#    print(ineq.tau)

#########
