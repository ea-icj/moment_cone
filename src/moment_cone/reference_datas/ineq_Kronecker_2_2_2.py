# Inequalities selected for V of KroneckerRepresentation type with dimensions [2, 2, 2, 1]
from moment_cone import *
G = LinearGroup([2, 2, 2, 1])
V = KroneckerRepresentation(G )
 
brut_inequations=[(1, 0, 1, 0, 0, 1, -2), 
 ] 

# inequalities in our formated type Inequality 
inequalities = [Inequality.from_tau(Tau.from_flatten(brut_ineq,G)) for brut_ineq in brut_inequations]
source = 'moment_cone'
