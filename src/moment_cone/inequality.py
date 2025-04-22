__all__ = (
    "Inequality",
    "full_under_symmetry_list_of_ineq",
)

from functools import cached_property
import itertools

from .typing import *
from .tau import *
from .permutation import Permutation
from .blocks import Blocks
from .root import Root
from .representation import *
from .rings import QQ, vector, Vector

class Inequality:
    """
    An inequality composed of a Tau and a list of permutations, one from each block of Tau.
    The list of permutation can be given as a list of permutations (w), as a list of inversions (inversions) 
    or as a dictionnary weight -> list of inversions. 

    In wtau, the blocks are permuted by the **inverse** of the corresponding permutation.
    
    Example :

    >>> from moment_cone import *
    >>> G = LinearGroup((4, 3, 2,1))
    >>> tau = Tau.from_flatten([6,2,1,4,1,2,5,3,1,1], G)
    >>> w = Permutation((1, 0, 3, 2)), Permutation((0, 2, 1)), Permutation((0, 1)),Permutation((0,))
    >>> ineq = Inequality(tau, w=w)
    >>> ineq
    Inequality(tau  = 6 2 1 4 | 1 2 5 | 3 1 | 1,
               w    = 1 0 3 2 | 0 2 1 | 0 1 | 0,
               wtau = 2 6 4 1 | 1 5 2 | 3 1 | 1)
    """
    tau: Tau

    def __init__(
        self, 
        tau: Tau, 
        w: Optional[Sequence[Permutation]] = None,
        inversions: Optional[Sequence[Root]] = None,
        gr_inversions: Optional[dict[int, list[Root]]] = None
    ):
        if w is None and inversions is None and gr_inversions is None :
            raise ValueError("At least one of w, inversions or gr_inversions, should be defined")
        assert w is None or len(tau.G) == len(list(w))
            
        self.tau = tau
        if w is not None:
            self.w = tuple(w)
        if inversions is not None:
            self.inversions = tuple(inversions)
        if gr_inversions is not None:
            self.gr_inversions = gr_inversions
    
    
    @cached_property
    def wtau(self) -> Tau:
        return Tau(tuple(wk.inverse(ck) for wk, ck in zip(self.w, self.tau.components)))
    
    @cached_property
    def w(self) -> tuple[Permutation, ...]:
        """
        Property w computed from inversions if not given at the initialisation.
        """
        return self.__compute_w()

    @cached_property
    def inversions(self) -> tuple[Root, ...]:
        """
        Property inversions computed from gr_inversions or w if not given at the initialisation.
        """
        if "gr_inversions" in self.__dict__:
            # Compute from gr_inversions
            return tuple(itertools.chain.from_iterable(self.gr_inversions.values()))
        else:
            # Compute from w
            return self.__compute_inversions_from_w()
        
    @cached_property
    def gr_inversions(self) -> dict[int, list[Root]]:
        """
        Property gr_inversions computed from inversions if not given at the initialisation.
        """
        return self.tau.grading_roots_in(self.inversions) 

    def __compute_w(self) -> tuple[Permutation, ...]:
        """
        Compute w from the inversions.
        """
        w=[]
        dims=self.tau.G
        inversions=self.inversions
        classified_inversions: list[list[Root]] = [[] for _ in range(len(dims))]
        for inv in inversions:
            classified_inversions[inv.k].append(inv)
        #print(classified_inversions)
        for l,d in enumerate(dims):
            i_inversions=[(inv.i,inv.j) for inv in classified_inversions[l]]
            w.append(Permutation.from_inversions(d,i_inversions))
        return tuple(w)

    def __compute_inversions_from_w(self) -> tuple[Root, ...]:
        """
        Compute the inversions from w.
        
        >>> from moment_cone import *
        >>> G = LinearGroup((2, 2, 2, 3, 1))
        >>> tau = Tau.from_flatten([6, 2, 1, 4, 1, 4, 5, 3, 1, 1], G)
        >>> w = Permutation((0, 1)), Permutation((1, 0)), Permutation((0, 1)), Permutation((2, 0, 1)),Permutation((0,))
        >>> ineq = Inequality(tau, w=w)
        >>> for r in ineq.inversions:
        ...     print(r)
        Root(k=1, i=0, j=1)
        Root(k=3, i=0, j=1)
        Root(k=3, i=0, j=2)
        """
        result=[]
        for k, p in enumerate(self.w):
            result += [Root(k,i,j) for i,j in p.inversions]
        return tuple(result)
        
    

    @staticmethod
    def from_tau(tau: Tau) -> "Inequality":
        """ Converts a (possibly non-dominant) tau to an element of the class Inequality

        An inequality is a pair (taup, w) where w.taup = tau and w is of minimal length with this property.
        
        >>> tau0 = Tau([[4, 9, 6, 5], [3, 1, 1, 2], [2, 2, 8, 2],[7]])
        >>> ineq0 = Inequality.from_tau(tau0)
        >>> Inequality.from_tau(tau0)
        Inequality(tau  = 9 6 5 4 | 3 2 1 1 | 8 2 2 2 | 7,
                   w    = 1 2 3 0 | 0 3 1 2 | 2 0 1 3 | 0,
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
        w = tuple((Permutation([i for t, i in taub]) for taub in tau_pairs))
        return Inequality(taup, w=w)
    
    @staticmethod
    def dominance(V: Representation, symmetry: bool = False) -> Iterable["Inequality"]:
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
             "           w    = " + " | ".join(" ".join(map(str, wk)) for wk in self.w) + ",\n" + \
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

        >>> from moment_cone import *
        >>> G = LinearGroup((2, 2, 2, 3, 1))
        >>> tau = Tau.from_flatten([6, 2, 1, 4, 1, 4, 5, 3, 1, 1], G)
        >>> w = Permutation((0, 1)), Permutation((1, 0)), Permutation((0, 1)), Permutation((2, 0, 1)), Permutation((0,))
        >>> ineq = Inequality(tau, w=w)
        >>> ineq
        Inequality(tau  = 6 2 | 1 4 | 1 4 | 5 3 1 | 1,
                   w    = 0 1 | 1 0 | 0 1 | 2 0 1 | 0,
                   wtau = 6 2 | 4 1 | 1 4 | 3 1 5 | 1)
        >>> ineq.sort_mod_sym_dim
        Inequality(tau  = 1 4 | 1 4 | 6 2 | 5 3 1 | 1,
                   w    = 0 1 | 1 0 | 0 1 | 2 0 1 | 0,
                   wtau = 1 4 | 4 1 | 6 2 | 3 1 5 | 1)
        """
        pairs = tuple(zip(self.tau.components, self.w))
        blocks = (sorted(b) for b in Blocks(pairs, self.tau.G.outer))
        tau_components, w = zip(*itertools.chain.from_iterable(blocks))
        tau = Tau(tau_components)
        return Inequality(tau, w=w)
    
    def weight_det(self,V: Representation) -> Vector:
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

