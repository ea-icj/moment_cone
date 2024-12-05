from .typing import *
from .tau import Tau
from .permutation import Permutation
from .blocks import Blocks

from functools import cached_property
import itertools

__all__ = (
    "Inequality",
)

class Inequality:
    """
    An inequality composed of a Tau and a permutation
    
    Example :
    >>> from cone import *
    >>> d = Dimension((2, 3, 4))
    >>> tau = Tau.from_flatten([1,6,2,1,4,1,2,5,3,1], d)
    >>> w = Permutation((1, 0)), Permutation((0, 2, 1)), Permutation((2, 0, 1, 3))
    >>> ineq = Inequality(tau, w)
    >>> print(ineq)    
    Inequality(tau  = 1 | 6 2 | 1 4 1 | 2 5 3 1,
               w    =     1 0 | 0 2 1 | 2 0 1 3,
               wtau = 1 | 2 6 | 1 1 4 | 3 2 5 1)
    """
    tau: Tau
    w: tuple[Permutation, ...]
    wtau: Tau # each w_k applied to each column of tau

    def __init__(self, tau: Tau, w: Iterable[Permutation]):
        self.tau = tau
        self.w = tuple(w)
        assert len(tau.d) == len(self.w)
        self.wtau = Tau(tuple(wk(ck) for wk, ck in zip(self.w, tau.components)), tau.ccomponent)

    def __repr__(self) -> str:
        return \
            f"Inequality(tau  = {self.tau},\n" + \
             "           w    =     " + " | ".join(" ".join(map(str, wk)) for wk in self.w) + ",\n" + \
            f"           wtau = {self.wtau})"
    
    @cached_property
    def sort_mod_sym_dim(self) -> "Inequality":
        """
        Sort (tau_i, w_i)_i by block of the dimensions

        >>> from cone import *
        >>> d = Dimension((2, 2, 2, 3))
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
        >>> w = Permutation((0, 1)), Permutation((1, 0)), Permutation((0, 1)), Permutation((2, 0, 1))
        >>> ineq = Inequality(tau, w)
        >>> ineq
        Inequality(tau  = 1 | 6 2 | 1 4 | 1 4 | 5 3 1,
                   w    =     0 1 | 1 0 | 0 1 | 2 0 1,
                   wtau = 1 | 6 2 | 4 1 | 1 4 | 1 5 3)
        >>> ineq.sort_mod_sym_dim
        Inequality(tau  = 1 | 1 4 | 1 4 | 6 2 | 5 3 1,
                   w    =     0 1 | 1 0 | 0 1 | 2 0 1,
                   wtau = 1 | 1 4 | 4 1 | 6 2 | 1 5 3)
        """
        pairs = tuple(zip(self.tau.components, self.w))
        blocks = (sorted(b) for b in Blocks(pairs, self.tau.d.symmetries))
        tau_components, w = zip(*itertools.chain.from_iterable(blocks))
        tau = Tau(tau_components, self.tau.ccomponent)
        return Inequality(tau, w)

