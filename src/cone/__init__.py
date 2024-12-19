#from sage.all import SymmetricFunctions, QQ, Integer as SageInteger
#sym_f = SymmetricFunctions(QQ).s()
#import sage.libs.lrcalc.lrcalc as lr

from .utils import *
from .blocks import *
from .partition import *
from .dimension import *
from .cone_dimension import *
from .tau import *
from .weight import *
from .root import *
from .task import *
from .inequality import *
from .permutation import *
from .hyperplane_candidates import *
from .rings import *
from .listW_temp import *
from .ramification import *
from .Normaliz2python import *
from .stabK_new2 import *
    

#def fct_Weights_of_Nu(nu: PartialMatrix[Partition]) -> PartialMatrix[int]:
#    """ Constructing size_nu from nu """
#    size_nu = PartialMatrix[int](*nu.shape)
#
#    for j, col in enumerate(nu.columns):
#        size_nu.extend(j, (cast(int, sum(p)) for p in col))
#
#    return size_nu



##### Pour 444, il y a un cas où on rentre dans LR_multi avec L=[] (liste vide) et nu=[].
## Cela bugue mais je ne comprends pas si c'est normal d'avoir une telle entrée


