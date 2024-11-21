#from sage.all import SymmetricFunctions, QQ, Integer as SageInteger
#sym_f = SymmetricFunctions(QQ).s()
#import sage.libs.lrcalc.lrcalc as lr

from .utils import *
from .typing import *
from .matrix import PartialMatrix
from .partition import Partition
from .dimension import Dimension
from .tau import Tau, ReducedTau
from .weight import Weight
from .root import Root
from .task import Task
from .inequality import Inequality

    

def fct_Weights_of_Nu(nu: PartialMatrix[Partition]) -> PartialMatrix[int]:
    """ Constructing size_nu from nu """
    size_nu = PartialMatrix[int](*nu.shape)

    for j, col in enumerate(nu.columns):
        size_nu.extend(j, (cast(int, sum(p)) for p in col))

    return size_nu



##### Pour 444, il y a un cas où on rentre dans LR_multi avec L=[] (liste vide) et nu=[].
## Cela bugue mais je ne comprends pas si c'est normal d'avoir une telle entrée
