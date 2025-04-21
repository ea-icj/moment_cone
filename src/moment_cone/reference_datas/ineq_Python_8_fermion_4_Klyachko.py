#Inequalities selected for V of boson type with dimensions [8]with number of particules = 4
__all__ = (
    'inequalities_K48',
)


from ..linear_group import LinearGroup
from ..representation import *
from ..inequality import Inequality
from ..tau import Tau

G=LinearGroup([8]) 
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
inequalities_K48=[Inequality.from_tau(Tau.from_flatten(brut_ineq,G)) for brut_ineq in brut_inequations]
