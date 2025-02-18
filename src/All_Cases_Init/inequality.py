from sage.all import QQ,vector
from functools import cached_property

import itertools

from .typing import *
from .tau import *
from .permutation import Permutation
from .blocks import Blocks
from .root import Root
from .rep import *

__all__ = (
    "Inequality",
    "full_under_symmetry_list_of_ineq",
)

class Inequality:
    """
    An inequality composed of a Tau and a list of permutations, one from each block of Tau.

    In wtau, the blocks are permuted by the **inverse** of the corresponding permutation.
    
    Example :
    >>> G = LinGroup((4, 3, 2,1))
    >>> tau = Tau.from_flatten([6,2,1,4,1,2,5,3,1,1], G)
    >>> w = Permutation((1, 0, 3, 2)), Permutation((0, 2, 1)), Permutation((0, 1)),Permutation((0,))
    >>> ineq = Inequality(tau, w)
    >>> Inequality(tau,w)
    Inequality(tau  = 6 2 1 4 | 1 2 5 | 3 1 | 1,
               w    = 1 0 3 2 | 0 2 1 | 0 1 | 0,
               wtau = 2 6 4 1 | 1 5 2 | 3 1 | 1)
    """
    tau: Tau
    w: tuple[Permutation, ...]
    wtau: Tau # each w_k applied to each column of tau

    def __init__(self, tau: Tau, w: Iterable[Permutation]):
        self.tau = tau
        self.w = tuple(w)
        assert len(tau.G) == len(self.w)
        self.wtau = Tau(tuple(wk.inverse(ck) for wk, ck in zip(self.w, tau.components)))
    
    @staticmethod
    def from_tau(tau: Tau) -> "Inequality":
        """ Converts a (possibly non-dominant) tau to an element of the class Inequality,
        that is a pair (taup, w) where w.taup = tau and w is of minimal length with this property.
        
        Example:
        >>> tau0 = Tau([[4, 9, 6, 5], [3, 1, 1, 2], [2, 2, 8, 2],[7]])
        >>> ineq0 = Inequality.from_tau(tau0)
        >>> Inequality.from_tau(tau0)
        Inequality(tau  = 9 6 5 4 | 3 2 1 1 | 8 2 2 2 | 7,
           w    =     1 2 3 0 | 0 3 1 2 | 2 0 1 3 | 0,
           wtau = 4 9 6 5 | 3 1 1 2 | 2 2 8 2 | 7)
        """
        tau_pairs = [
            sorted(
                ((t, i) for i, t in enumerate(taub)),
                key=lambda pair: (-pair[0], pair[1])
            )
            for taub in tau.components
        ]

        taup = Tau(
            Blocks.from_blocks([[t for t, i in taub] for taub in tau_pairs])
        )
        w = (Permutation([i for t, i in taub]) for taub in tau_pairs)
        return Inequality(taup, w)
    
    @staticmethod
    def dominance(V: Representation, symmetry=False) -> Iterable["Inequality"]:
        """
        Computes the dominant inequalities 
        if symmetry=True, only inequalities up to symmetries of G are computed.
        """
        Res=[]
        for k,dk in enumerate(V.G):
            for i in range(dk-1):
                component=i*[0]+[-1,1]+(dk-i-2)*[0]
                tau=Tau([dj*[0] for dj in V.G[:k]]+[component]+[dj*[0] for dj in V.G[k+1:]])
                Res.append(Inequality.from_tau(tau))
            #if V.type=='kron' and k!=len(V.G)-1:
            #    component=(dk-1)*[0]+[-1]
            #    tau=Tau([dj*[0] for dj in V.G[:k]]+[component]+[dj*[0] for dj in V.G[k+1:]],V.G)
            #    Res.append(Inequality.from_tau(tau))
        if not(symmetry):
            return(Res)
        else:
            return(list(set(ineq.sort_mod_sym_dim for ineq in Res)))

    
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

        >>> from All_Init import *
        >>> G = LinGroup((2, 2, 2, 3))
        >>> tau = Tau.from_flatten([6, 2, 1, 4, 1, 4, 5, 3, 1, 1], G)
        >>> w = Permutation((0, 1)), Permutation((1, 0)), Permutation((0, 1)), Permutation((2, 0, 1)), Permutation((0,))
        >>> ineq = Inequality(tau, w)
        >>> ineq
        Inequality(tau  = 1 | 6 2 | 1 4 | 1 4 | 5 3 1,
                   w    =     0 1 | 1 0 | 0 1 | 2 0 1,
                   wtau = 1 | 6 2 | 4 1 | 1 4 | 3 1 5)
        >>> ineq.sort_mod_sym_dim
        Inequality(tau  = 1 | 1 4 | 1 4 | 6 2 | 5 3 1,
                   w    =     0 1 | 1 0 | 0 1 | 2 0 1,
                   wtau = 1 | 1 4 | 4 1 | 6 2 | 3 1 5)
        """
        pairs = tuple(zip(self.tau.components, self.w))
        blocks = (sorted(b) for b in Blocks(pairs, self.tau.G.outer))
        tau_components, w = zip(*itertools.chain.from_iterable(blocks))
        tau = Tau(tau_components)
        return Inequality(tau, w)

    @property
    def inversions(self) -> Iterable[Root]:
        """
        Returns all possible inversions Root(k, i, j) of w
        
        >>> from cone import *
        >>> G = LinGroup((2, 2, 2, 3,1))
        >>> tau = Tau.from_flatten([6, 2, 1, 4, 1, 4, 5, 3, 1,1], G)
        >>> w = Permutation((0, 1)), Permutation((1, 0)), Permutation((0, 1)), Permutation((2, 0, 1)),Permutation((0,))
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


    
    def weight_det(self,V: Representation) -> vector:
        """
        Weight chi_det of Theorem BKR
        """
        tau=self.tau
        listp=[]
        for ll in list(tau.positive_weights(V).values()):
            listp+=ll
        if listp == [] and list(self.inversions)==[]:
            return(vector(QQ,sum(V.G)))
        else :
            return(sum([chi.as_vector for chi in listp])-sum([root.to_vector(V.G) for root in self.inversions]))
            

def full_under_symmetry_list_of_ineq(seq_ineq: Iterable[Inequality]) -> Iterable[Inequality] :
    seq_tau=full_under_symmetry_list_of_tau([ineq.wtau for ineq in seq_ineq])
    return([Inequality.from_tau(tau) for tau in seq_tau])


