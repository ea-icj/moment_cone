#Inequalities selected for V of fermion type with dimensions [7]with number of particules = 3
from moment_cone import *

G = LinearGroup([7]) 
V = FermionRepresentation(G, particle_cnt=3 )
 
brut_inequations=[(-2, 1, 1, 1, 1, -2, -2), 
(1, -2, 1, 1, -2, 1, -2), 
(1, 1, -2, -2, 1, 1, -2), 
(1, 1, -2, 1, -2, -2, 1), 
 ] 

#inequalities in our formated type Inequality 
inequalities=[Inequality.from_tau(Tau.from_flatten(brut_ineq,G)) for brut_ineq in brut_inequations] 
source = "Klyachko"
 
