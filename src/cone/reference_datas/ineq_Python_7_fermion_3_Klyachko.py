#Inequalities selected for V of fermion type with dimensions [7]with number of particules = 3
__all__ = (
    'inequalities_K37',
)

from ..linear_group import LinearGroup
from ..representation import *
from ..inequality import Inequality
from ..tau import Tau

G=LinearGroup([7]) 
V = FermionRepresentation(G, particle_cnt=3 )
 
brut_inequations=[(-2, 1, 1, 1, 1, -2, -2), 
(1, -2, 1, 1, -2, 1, -2), 
(1, 1, -2, -2, 1, 1, -2), 
(1, 1, -2, 1, -2, -2, 1), 
 ] 

#inequalities in our formated type Inequality 
inequalities_K37=[Inequality.from_tau(Tau.from_flatten(brut_ineq,G)) for brut_ineq in brut_inequations] 
 
 
