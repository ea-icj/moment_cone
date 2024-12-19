from .typing import *
from .dimension import Dimension
from .partition import *
from .blocks import Blocks
from .weight import Weight
from .root import Root
from .rings import matrix, ZZ, QQ
from .utils import extend_with_repetitions, flatten_dictionary, grading_dictionary, orbit_symmetries
from .permutation import Permutation

import itertools
from functools import cached_property

"""
TODO:
- same method for both scalar product using @singledispatch?
- property or method for Pzero?
"""

__all__ = (
    "Tau",
    "ReducedTau",
    "unique_modulo_symmetry_list_of_tau",
    "full_under_symmetry_list_of_tau",
    "find_1PS_reg_mod_sym_dim",
    "find_1PS_mod_sym_dim",
)

class Tau:
    """
    Tuple of partition along with a coefficient
    
    Example:
    >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)), 1)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1

    >>> d = Dimension((4, 3, 3))
    >>> tau = Tau.from_flatten((1, 3, 3, 2, 2, 2, 2, 1, 2, 2, 1), d)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1
    """
    #__slots__ = 'ccomponent', '_components' # FIXME cached_property cannot work without __dict__ ... => self managed cache or removing __slots__
    ccomponent: Optional[int]
    _components: Blocks[int]

    def __init__(self, components: Iterable[Sequence[int]] | Blocks[int], ccomponent: Optional[int] = None):
        """
        Tau initialization from a sequence of sub-group or directly from a partial matrix
        
        Beware that the ccomponent is at last position in the arguments since it is optional.
        """
        self.ccomponent = ccomponent
        if isinstance(components, Blocks):
            self._components = components.freeze()
        else:
            self._components = Blocks.from_blocks(components)

    @staticmethod
    def from_flatten(s: Iterable[int], d: Dimension) -> "Tau":
        """ Returns tau from its flattened version """
        all_components = tuple(s)
        if len(all_components) == d.sum:
            ccomponent = None
        elif len(all_components) == d.sum + 1:
            ccomponent = all_components[0]
            all_components = all_components[1:]
        else:
            raise ValueError("Invalid number of components")

        # Tau will be always immutable
        return Tau(Blocks.from_flatten(tuple(all_components), d), ccomponent)

    @staticmethod
    def from_zero_weights(weights: Sequence[Weight], d: Dimension) -> "Tau":
        """
        From a set of weights generating an hyperplane in X^*(T), returns a primitive Tau orthogonal to the hyperplane

        TODO: doctest
        """
        from .hyperplane_candidates import hyperplane_matrix

        M = hyperplane_matrix(weights, d)
        Mp = M.augment(matrix(ZZ, [len(d) * [0] for i in range(d.sum + 1)])) 
        for u in range(len(d)):
           # FIXME: typically the case of usage of Blocks but it may not work with Sage matrix ...
           shift = sum(d[k] for k in range(u))
           for i in range(d[u]):
               Mp[shift + i + 1,u + len(weights)] = 1

        b = Mp.kernel().basis()
        if len(b) != 1:
           raise ValueError("Given set of weights does not generates an hyperplane")
        else:
           return Tau.from_flatten(b[0], d)

    @property
    def opposite(self) -> "Tau":
        """
        Returns the opposite of some tau (same orthogonal hyperplane)

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> tau.opposite
        7 | -3 -2 -2 | -4 -2 -1 | -3 -2
        """
        return Tau.from_flatten(
            (-x for x in self.flattened),
            self.d
        )

    def __len__(self) -> int:
        """ Number of components """
        return len(self._components)

    @cached_property
    def d(self) -> Dimension:
        """ Length of each component """
        return Dimension(self._components.sizes)

    @cached_property
    def components(self) -> tuple[tuple[int, ...], ...]:
        """ Sequence of the components of tau """
        return cast(tuple[tuple[int, ...], ...], tuple(self._components.blocks))

    @property
    def flattened(self) -> Iterable[int]:
        """
        Returns the whole tau as a unique sequence.

        ie (cc-component | column1 | column2 | ...)
        """
        if self.ccomponent is None:
            return self._components.flatten
        else:
            return (self.ccomponent,) + cast(tuple[int], self._components.flatten)

    @cached_property
    def reduced(self) -> "ReducedTau":
        """ Returns reduced form of tau """
        return ReducedTau(self)
        
    def m_extend_with_repetitions(self,d:Dimension)-> Iterable["Tau"]:
        """
        Extends self in all possible manners to a Tau of dimension d
                
        Example:
        >>> tau = Tau([[5, 4, 0], [4, 3, 0], [4, 2, 0]], -1)
        >>> tau
        -1 | 5 4 0 | 4 3 0 | 4 2 0
        >>> for t in tau.m_extend_with_repetitions(Dimension([4, 4, 3])):
        ...     print(t)
        -1 | 5 4 0 0 | 4 3 0 0 | 4 2 0
        -1 | 5 4 0 0 | 4 3 3 0 | 4 2 0
        -1 | 5 4 0 0 | 4 4 3 0 | 4 2 0
        -1 | 5 4 4 0 | 4 3 0 0 | 4 2 0
        -1 | 5 4 4 0 | 4 3 3 0 | 4 2 0
        -1 | 5 4 4 0 | 4 4 3 0 | 4 2 0
        -1 | 5 5 4 0 | 4 3 0 0 | 4 2 0
        -1 | 5 5 4 0 | 4 3 3 0 | 4 2 0
        -1 | 5 5 4 0 | 4 4 3 0 | 4 2 0
        """
        extend_each_comp=[extend_with_repetitions(x,d[i]) for i,x in enumerate(self._components)]
        for ext in itertools.product(*extend_each_comp):
            yield Tau(ext,self.ccomponent)

    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + repr(self._components)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Tau):
            return NotImplemented
        return self.ccomponent == other.ccomponent and self._components == other._components

    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash((self.ccomponent, self._components))

    # TODO: or same dot method with @singledispatch ?
    def dot_weight(self, weight: Weight, ccomponent: Optional[int] = None) -> int:
        """ Scalar product of tau with a weight of V """
        assert not (self.ccomponent is None and ccomponent is None)
        if ccomponent is None:
            ccomponent = self.ccomponent
        return cast(int, ccomponent) + sum(c[wi] for c, wi in zip(self.components, weight))

    def dot_root(self, root: Root) -> int:
        """ Scalar product of tau with a root of U """
        c = self.components[root.k]
        return c[root.i] - c[root.j]

    @cached_property
    def is_dom_reg(self) -> bool:
        """ Check if tau is dominant and regular """
        return all(all(a > b for a, b in itertools.pairwise(c)) for c in self.components)

    @cached_property
    def is_dominant(self) -> bool:
        """ Check if tau is dominant """
        return all(all(a >= b for a, b in itertools.pairwise(c)) for c in self.components)
    
    @cached_property
    def sl_representative(self) -> "Tau":
        """ Returns representative of tau in C^* x (SLn)^3 """
        assert self.ccomponent is not None
        from math import lcm, gcd
        #flat_non_zero=[x for x in self.flattened if x!=0]
        tau_lcm = lcm(*self.d)
        ccomponent = self.ccomponent * tau_lcm
        columns = []
        for dj, cj in zip(self.d, self.components):
            column_sum = sum(cj)
            shift = column_sum * tau_lcm // dj
            columns.append([tau_lcm * cji - shift for cji in cj])
            ccomponent += shift
        res_gcd = gcd(ccomponent, *itertools.chain.from_iterable(columns))
        return Tau(
            tuple(tuple(v // res_gcd for v in cj) for cj in columns),
            ccomponent // res_gcd
        )

    @cached_property
    def end0_representative(self) -> "Tau":
        """ Returns representative of tau in ??? space (final value of each block is zero) """
        assert self.ccomponent is not None
        from math import gcd
        total_shift = 0
        columns = []
        for cj in self.components:
            columns.append([cji - cj[-1] for cji in cj])
            total_shift += cj[-1]
        ccomponent = self.ccomponent + total_shift

        res_gcd = gcd(ccomponent, *itertools.chain.from_iterable(columns))
        return Tau(
            tuple(tuple(v // res_gcd for v in cj) for cj in columns),
            ccomponent // res_gcd
        )

    def grading_weights_in(self, weights: Iterable[Weight]) -> dict[int, list[Weight]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on a give subset of V.
        """
        from .utils import grading_dictionary
        return grading_dictionary(weights, self.dot_weight)

    @cached_property
    def grading_weights(self) -> dict[int, list[Weight]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on V.
        
        For each key p, the weights in the entry p correspond to a basis of the eigen space
        
        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.grading_weights
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        -2: [Weight((1, 2, 1), idx: 11), Weight((2, 2, 1), idx: 17)]
        -1: [Weight((0, 2, 1), idx: 5), Weight((1, 1, 1), idx: 9), Weight((1, 2, 0), idx: 10), Weight((2, 1, 1), idx: 15), Weight((2, 2, 0), idx: 16)]
        0: [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]
        1: [Weight((0, 1, 0), idx: 2), Weight((1, 0, 1), idx: 7), Weight((2, 0, 1), idx: 13)]
        2: [Weight((0, 0, 1), idx: 1), Weight((1, 0, 0), idx: 6), Weight((2, 0, 0), idx: 12)]
        3: [Weight((0, 0, 0), idx: 0)]
        """
        return self.grading_weights_in(Weight.all(self.d))

    def grading_roots_in(self, roots: Iterable[Root]) -> dict[int, list[Root]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on the given subset of U.
        """
        from .utils import grading_dictionary
        return grading_dictionary(roots, self.dot_root)


    @cached_property
    def grading_roots(self) -> dict[int, list[Root]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on U (sum of positive root spaces).
        
        For each key p, the roots in the entry p correspond to a basis of the eigen space.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gr = tau.grading_roots
        >>> for k in sorted(gr.keys()):
        ...     print(f"{k}:", gr[k])
        0: [Root(k=0, i=1, j=2)]
        1: [Root(k=0, i=0, j=1), Root(k=0, i=0, j=2), Root(k=1, i=1, j=2), Root(k=2, i=0, j=1)]
        2: [Root(k=1, i=0, j=1)]
        3: [Root(k=1, i=0, j=2)]
        """
        return self.grading_roots_in(Root.all_of_U(self.d))

    @property
    def positive_weights(self) -> dict[int, list[Weight]]:
        """
        Basis of the eigen space for positive eigen values for the action of tau on V.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.positive_weights
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        1: [Weight((0, 1, 0), idx: 2), Weight((1, 0, 1), idx: 7), Weight((2, 0, 1), idx: 13)]
        2: [Weight((0, 0, 1), idx: 1), Weight((1, 0, 0), idx: 6), Weight((2, 0, 0), idx: 12)]
        3: [Weight((0, 0, 0), idx: 0)]
        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_weights, lambda x: x > 0)

    @property
    def non_positive_weights(self) -> dict[int, list[Weight]]:
        """
        Basis of the eigen space for non-positive eigen values for the action of tau on V.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.non_positive_weights
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        -2: [Weight((1, 2, 1), idx: 11), Weight((2, 2, 1), idx: 17)]
        -1: [Weight((0, 2, 1), idx: 5), Weight((1, 1, 1), idx: 9), Weight((1, 2, 0), idx: 10), Weight((2, 1, 1), idx: 15), Weight((2, 2, 0), idx: 16)]
        0: [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]

        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_weights, lambda x: x <= 0)


    @property
    def non_negative_weights(self) -> dict[int, list[Weight]]:
        """
        Basis of the eigen space for non-negative eigen values for the action of tau on V.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.non_negative_weights
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        0: [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]
        1: [Weight((0, 1, 0), idx: 2), Weight((1, 0, 1), idx: 7), Weight((2, 0, 1), idx: 13)]
        2: [Weight((0, 0, 1), idx: 1), Weight((1, 0, 0), idx: 6), Weight((2, 0, 0), idx: 12)]
        3: [Weight((0, 0, 0), idx: 0)]
        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_weights, lambda x: x >= 0)
    
    @property
    def positive_roots(self) -> dict[int, list[Root]]:
        """
        Basis of the eigen space for positive eigen values for the action of tau on U.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gr = tau.positive_roots
        >>> for k in sorted(gr.keys()):
        ...     print(f"{k}:", gr[k])
        1: [Root(k=0, i=0, j=1), Root(k=0, i=0, j=2), Root(k=1, i=1, j=2), Root(k=2, i=0, j=1)]
        2: [Root(k=1, i=0, j=1)]
        3: [Root(k=1, i=0, j=2)]
        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_roots, lambda x: x > 0)

    #TODO: definition changed thanks to above definition.
    #Optimization might be done to avoid lists and dict, dealing only with iterables?
    @property
    def orthogonal_roots(self) -> list[Root]:
        """
        All the positive roots beta of V so that <beta, tau> = 0

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> tau.orthogonal_roots
        [Root(k=0, i=1, j=2)]
        """
        return self.grading_roots.get(0, [])

    @property
    def orthogonal_all_roots(self) -> list[Root]:
        """
        All the roots beta of V so that <beta, tau> = 0. beta=Root(k,i,i) allowed. 
        We get a set indexing a bases of Lie(K^tau)

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> for t in tau.orthogonal_all_roots:
        ...     print(t)
        Root(k=0, i=0, j=0)
        Root(k=0, i=1, j=1)
        Root(k=0, i=1, j=2)
        Root(k=0, i=2, j=1)
        Root(k=0, i=2, j=2)
        Root(k=1, i=0, j=0)
        Root(k=1, i=1, j=1)
        Root(k=1, i=2, j=2)
        Root(k=2, i=0, j=0)
        Root(k=2, i=1, j=1)
        """
        res=[]
        for alpha in Root.all_of_K(self.d):
            if self.dot_root(alpha)==0:
                res.append(alpha)
        return(res)        

    @property
    def orthogonal_weights(self) -> list[Weight]:
        """
        All the root beta so that <beta, tau> = 0
        
        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> tau.orthogonal_weights
        [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]
        """
        return self.grading_weights.get(0, [])

    @cached_property
    def sort_mod_sym_dim(self) -> "Tau":
        """
        Sort tau by block of the dimensions
        
        >>> d = Dimension((2, 2, 2, 1, 1, 1))
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
        >>> tau
        1 | 6 2 | 1 4 | 1 4 | 5 | 3 | 1
        >>> tau.sort_mod_sym_dim
        1 | 1 4 | 1 4 | 6 2 | 1 | 3 | 5
        """
        blocks = (sorted(b) for b in Blocks(self.components, self.d.symmetries))
        return Tau(itertools.chain.from_iterable(blocks), self.ccomponent)
        
    def orbit_symmetries(self) -> Iterable["Tau"]:
        """
        Lists the orbit of tau under symmetries of dimensions of its components

        Example:
        >>> d = Dimension((2, 2, 2, 1, 1))
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3], d)
        >>> tau
        1 | 6 2 | 1 4 | 1 4 | 5 | 3
        >>> for t in tau.orbit_symmetries():
        ...     print(t)
        1 | 1 4 | 1 4 | 6 2 | 3 | 5
        1 | 1 4 | 1 4 | 6 2 | 5 | 3
        1 | 1 4 | 6 2 | 1 4 | 3 | 5
        1 | 1 4 | 6 2 | 1 4 | 5 | 3
        1 | 6 2 | 1 4 | 1 4 | 3 | 5
        1 | 6 2 | 1 4 | 1 4 | 5 | 3
        """
        for sym_comp in orbit_symmetries(self._components, self.d.symmetries):
            yield Tau(sym_comp, self.ccomponent)

    @cached_property
    def is_sub_module(self) -> bool :
        """
        Check if V^{\tau>0} is a C^*-submodule of U

        TODO: doctest
        """
        # Create the dictionnary for V^{\tau>0}
        gw = self.positive_weights 
        # Create the dictionnary for U
        # FIXME: it's self.grading_roots
        gr = grading_dictionary(Root.all_of_U(self.d), self.dot_root)

        # Compare them
        # TODO: using all instead (after doctest) ?
        for p in gw.keys():
            if p not in gr.keys() or len(gr[p])<len(gw[p]) :
                return False
        return True

    @cached_property
    def dim_Pu(self) -> int:
        """
        Dimension of Pu
        
        Example:
        >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)), 1)
        >>> tau
        1 | 3 3 2 2 | 2 2 1 | 2 2 1
        >>> tau.dim_Pu
        8
        >>> sum(len(roots) for roots in tau.positive_roots.values())
        8
        """
        sum_di2 = sum(di**2 for di in self.d)
        sum_mi2 = sum(mi**2 for mi in self.reduced.mult.flatten)
        return (sum_di2 - sum_mi2) // 2

class ReducedTau:
    """
    Tau in a reduction form
    
    Example:
    >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)), 1)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1
    >>> rtau = tau.reduced
    >>> rtau
    1 | 3^2 2^2 | 2^2 1^1 | 2^2 1^1
    >>> rtau.values
    3 2 | 2 1 | 2 1
    >>> rtau.mult
    2 2 | 2 1 | 2 1
    """
    __slots__ = 'ccomponent', 'values', 'mult'
    ccomponent: Optional[int]
    values: Blocks[int]
    mult: Blocks[int]

    def __init__(self, tau: Tau):
        from .utils import group_by_block
        values = Blocks.from_flatten([], (0,) * len(tau.d))
        mult = Blocks.from_flatten([], (0,) * len(tau.d))

        for j, component in enumerate(tau.components):
            values[j], mult[j] = zip(*group_by_block(component))

        self.ccomponent = tau.ccomponent
        self.values = values.freeze()
        self.mult = mult.freeze()

    def __len__(self) -> int:
        """ Number of components """
        return len(self.values)

    @property
    def small_d(self) -> Dimension:
        return Dimension(self.values.sizes)

    # FIXME: it seems that i,j is not really needed => simplification?
    def __getitem__(self, idx: tuple[int, int]) -> tuple[int, int]:
        return self.values[idx[1], idx[0]], self.mult[idx[1], idx[0]]

    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + " | ".join(
            " ".join(f"{v}^{m}" for v, m in zip(cv, cm))
            for cv, cm in zip(self.values.blocks, self.mult.blocks)
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ReducedTau):
            return NotImplemented
        return self.ccomponent == other.ccomponent and self.values == other.values and self.mult == other.mult

    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash((self.ccomponent, self.values, self.mult))

    # FIXME: returns an generator and not a sequence. Is it adapted to a property (or a cached_property?)
    # TODO: it will be used many time => returns a sequence and @cached_property
    # TODO: remove? information already contained in orthogonal weight like subsequent function.
    @property
    def Pzero(self) -> Iterable[Weight]:
        """ Search for weights w of V so that C_component + sum_k tau_red[w_k, k] = 0 """
        as_tau = Tau(self.values, self.ccomponent) # We use the scalar product of Tau
        for weight in Weight.all(self.small_d):
            if as_tau.dot_weight(weight) == 0:
                yield weight

def unique_modulo_symmetry_list_of_tau(seq_tau: Iterable[Tau]) -> set[Tau]:
    """
    Unique sequence of tau modulo the it's symmetries

    Example:
    >>> d = Dimension((2, 2, 2, 1, 1))
    >>> t1 = Tau.from_flatten([1, 6, 2, 1, 5, 1, 4, 5, 3], d)
    >>> t2 = t1.end0_representative
    >>> t3 = t2.sort_mod_sym_dim
    >>> t4 = t1.sort_mod_sym_dim
    >>> t5 = Tau.from_flatten([2, 6, 2, 1, 5, 1, 4, 5, 3], d)
    >>> for tau in unique_modulo_symmetry_list_of_tau((t1, t2, t3, t4, t5)):
    ...     print(tau)
    21 | -4 0 | -3 0 | 4 0 | 0 | 0
    20 | -4 0 | -3 0 | 4 0 | 0 | 0
    """
    return {tau.end0_representative.sort_mod_sym_dim for tau in seq_tau}

def full_under_symmetry_list_of_tau(seq_tau: Iterable[Tau]) -> Iterable[Tau]:
    """ TODO """
    return itertools.chain.from_iterable(tau.orbit_symmetries() for tau in seq_tau)

def find_1PS_reg_mod_sym_dim(d: Dimension, u: int) -> Iterable[Tau]:
    """
    Lists candidates dominant regular 1 parameter subgroups tau, whose orthogonal is generated by weights, and such that there is at most u weights we of V with tau(we)>0
    With results up to the action of the symmetries of d.

    TODO: move to a dedicated file?

    Example:
    >>> d = Dimension((4, 4, 3))
    >>> for tau in find_1PS_reg_mod_sym_dim(d, 15):
    ...     print(tau)
    -5 | 3 2 1 0 | 3 2 1 0 | 3 1 0
    -6 | 3 2 1 0 | 3 2 1 0 | 3 2 0
    -7 | 4 3 2 0 | 4 3 2 0 | 3 1 0
    -6 | 4 3 2 0 | 3 2 1 0 | 3 1 0
    -5 | 4 2 1 0 | 3 2 1 0 | 2 1 0
    -5 | 3 2 1 0 | 3 2 1 0 | 2 1 0
    -6 | 4 3 2 0 | 4 2 1 0 | 2 1 0
    """
    from .hyperplane_candidates import find_hyperplanes_reg_mod_sym_dim
    list_hr = find_hyperplanes_reg_mod_sym_dim(d, u)
    list_1PS = unique_modulo_symmetry_list_of_tau(Tau.from_zero_weights(h, d) for h in list_hr)
    list_1PS_sign = (t for tau in list_1PS for t in (tau, tau.opposite))
    return filter(lambda tau: tau.is_dom_reg, list_1PS_sign)


def find_1PS_mod_sym_dim(d: Dimension) -> Sequence["Tau"]:
    """
    Same as find_1PS_reg_mod_sym_dim without regularity condition
    Computed by 

    TODO: move to a dedicated file?    
    """
    # Initialisation with regular 1-PS
    Liste_1PS = list(find_1PS_reg_mod_sym_dim(d,d.u_max(d)))
    print('For d=',d,'we get',len(Liste_1PS),' candidates regular dominant up to symmetry')
    # Looking for 1-PS by extension
    sub_dim=[Dimension(p) for p in Partition(d).all_subpartitions()][1:-1] #[1:-1] excludes 1... 1 and d
    for small_d in sub_dim:  
        umax=d.u_max(small_d)
        #Recover by induction all candidates 1-PS mod symmetry
        Liste_1PS_smalld_reg_mod_sym= list(find_1PS_reg_mod_sym_dim(small_d,umax))
        print('For d=',small_d,'we get',len(Liste_1PS_smalld_reg_mod_sym),' candidates regular dominant up to symmetry')
        Liste_1PS_smalld_reg=full_under_symmetry_list_of_tau(Liste_1PS_smalld_reg_mod_sym)
        #Liste_1PS_smalld_reg=sum([list(tau.orbit_symmetries()) for tau in Liste_1PS_smalld_reg_mod_sym]  ,[])
        Liste_1PS_smalld_extended=[]
        for permut in Permutation.embeddings_mod_sym(d, small_d):
            for tau in Liste_1PS_smalld_reg:
                tau_twist=Tau([tau._components[i] for i in permut],tau.ccomponent)
                list_tau_extended=tau_twist.m_extend_with_repetitions(d)
                for tau_ext in list_tau_extended:
                    if len(flatten_dictionary(tau_ext.positive_weights))<=tau_ext.dim_Pu:
                        Liste_1PS_smalld_extended.append(tau_ext)
        Liste_1PS+=unique_modulo_symmetry_list_of_tau(Liste_1PS_smalld_extended)
    return(Liste_1PS)



#    return list(set(Liste_1PS))
