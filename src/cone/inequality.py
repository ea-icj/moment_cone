from .typing import *
from .tau import Tau
from .permutation import Permutation
from .blocks import Blocks
from .root import Root


from functools import cached_property
import itertools

__all__ = (
    "Inequality",
    "unique_modulo_symmetry_list_of_ineq",
)

class Inequality:
    """
    An inequality composed of a Tau and a list of permutations
    
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
   
    def from_tau(tau:Tau) -> "Inequality":
        """converts a (possibly non-dominant) tau to an element of the class Inequality, that is a pair (taup,w) where w.taup=tau
        
        """
        tau_pairs=[sorted([(t,i) for i,t in enumerate(taub)],key=lambda pair:(-pair[0],pair[1])) for taub in tau._components]
        for taub in tau_pairs:
            taup=Tau(Blocks.from_blocks([[t for t,i in taub] for taub in tau_pairs]), tau.ccomponent)
            w=[Permutation([i for t,i in taub]).inverse for taub in tau_pairs]
        return Inequality(taup,w)
    
    def __repr__(self) -> str:
        return \
            f"Inequality(tau  = {self.tau},\n" + \
             "           w    =     " + " | ".join(" ".join(map(str, wk)) for wk in self.w) + ",\n" + \
            f"           wtau = {self.wtau})"
    
    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Inequality):
            return NotImplemented
        return self.tau == other.tau and self.w == other.w
    
    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash((self.tau, self.w))
    
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

    @property
    def inversions(self) -> Iterable[Root]:
        """
        Returns all possible inversions Root(k, i, j) of w
        
        >>> from cone import *
        >>> d = Dimension((2, 2, 2, 3))
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
        >>> w = Permutation((0, 1)), Permutation((1, 0)), Permutation((0, 1)), Permutation((2, 0, 1))
        >>> ineq = Inequality(tau, w)
        >>> for r in ineq.inversions:
        ...     print(r)
        Root(k=1, i=0, j=1)
        Root(k=3, i=0, j=1)
        Root(k=3, i=0, j=2)
        """
        for k, p in enumerate(self.w):
            for i, j in p.inversions:
                yield Root(k, i, j)

def unique_modulo_symmetry_list_of_ineq(seq_ineq:Sequence["Inequality"]) -> Sequence["Inequality"]:
    return list(set([ineq.sort_mod_sym_dim for ineq in seq_ineq]))


