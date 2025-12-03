#Inequalities selected for V of boson type with dimensions [8]with number of particules = 4
from moment_cone import *

G = LinearGroup([8]) 
V = FermionRepresentation(G, particle_cnt=4 )
 
brut_inequations=[(1, 0, 0, -1, 0, -1, -1, 0), 
(0, 0, 1, -1, 0, 0, -1, -1), 
(0, 0, 0, 0, 1, -1, -1, -1), 
(0, 1, 0, -1, 0, -1, 0, -1), 
(1, -1, 0, 0, 0, 0, -1, -1), 
(1, 0, 0, -1, -1, 0, 0, -1), 
(1, 0, -1, 0, 0, -1, 0, -1), 
(-1, 1, 1, -1, 1, -1, -1, -3), 
(1, 1, -1, -1, 1, -3, -1, -1), 
(1, 1, 1, -3, -1, -1, -1, -1), 
(1, -1, 1, -1, 1, -1, -3, -1), 
(1, 1, -1, -1, -1, -1, 1, -3), 
(1, -1, -1, 1, 1, -1, -1, -3), 
(1, -1, 1, -1, -1, 1, -1, -3), 
(3, -1, -1, -1, -1, -1, -1, -1), 
 ] 

#inequalities in our formated type Inequality 
inequalities=[Inequality.from_tau(Tau.from_flatten(brut_ineq,G)) for brut_ineq in brut_inequations]
source = "Klyachko"
