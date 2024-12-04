from .typing import *
from .dimension import Dimension
from .blocks import Blocks
from .weight import Weight
from .root import Root
from .hyperplane_candidates import hyperplane_matrix
from sage.all import matrix, ZZ, QQ # type: ignore

import itertools
from functools import cached_property

"""
TODO:
- test file
- same method for both scalar product using @singledispatch?
- property or method for Pzero?
"""

__all__ = (
    "Tau",
    "ReducedTau"
)

class Tau:
    """ Tuple of partition along with a coefficient """
    #__slots__ = 'ccomponent', '_components' # FIXME cached_property cannot work without __dict__ ... => self managed cache or removing __slots__
    ccomponent: Optional[int]
    _components: Blocks[int]

    def __init__(self, components: Iterable[Sequence[int]] | Blocks[int], ccomponent: Optional[int] = None):
        """ Tau initialization from a sequence of sub-group or directly from a partial matrix """
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
    def from_zero_weights(S: Sequence[Weight], d: Dimension) -> "Tau":
        """ From a set of weights generating an hyperplane in X^*(T), returns a primitive Tau orthogonal to the hyperplane"""
        M = hyperplane_matrix(S, d)
        MQ = M.change_ring(QQ)
        MQ.augment(matrix(QQ,[len(d)*[0] for i in range(sum(d)+1)])) 
        for u in range(len(d)):
           shift=sum(d[k] for k in range(u))
           for i in range(d[u]):
               M[shift+i+1,u+len(S)]=1
        b=M.kernel().basis()
        if len(b)!=1:
           raise ValueError("Given set of weights does not generates an hyperplane")
        else:
           return Tau.from_flatten(b[0],d)

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
        tau_lcm = lcm(*self.flattened)
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

    # TODO: generate the dictionary for all values of the product scalar
    # and filtering it later. Renaming it like grading_weights and removing
    # the optional weights list so that to be a @cached_property.
    def positive_weights(self, weights: Optional[Iterable[Weight]] = None) -> dict[int, list[Weight]]:
        """ Inverse image of each non-negative p = <w, tau> for each w in weights (all weights by default) """
        if weights is None:
            weights = Weight.all(self.d)

        result: dict[int, list[Weight]] = {}
        for chi in weights:
            p = self.dot_weight(chi)
            if p >= 0:
                result.setdefault(p, []).append(chi)
        return result

    # TODO: generate the dictionary for all values of the product scalar
    # and filtering it later. Renaming it like grading_roots and removing
    # the optional weights list so that to be a @cached_property.
    def positive_roots(self, roots: Optional[Iterable[Root]] = None) -> dict[int, list[Root]]:
        """ Inverse image of each non-zero p = <beta, tau> for each beta in roots (all roots of U by default) """
        if roots is None:
            roots = Root.all(self.d)

        result: dict[int, list[Root]] = {}
        for r in roots:
            p = self.dot_root(r)
            if p > 0:
                result.setdefault(p, []).append(r)
        return result

    # TODO: as a @cache_property ?
    def orthogonal_roots(self) -> Iterable[Root]:
        """ All the root beta so that <beta, tau> = 0 """
        return filter(
            lambda root: self.dot_root(root) == 0,
            Root.all(self.d)
        )

    @cached_property
    def sort_mod_sym_dim(self) -> "Tau":
        """ Sort tau by block of the dimensions """
        blocks = (sorted(b) for b in Blocks(self.components, self.d.symmetries))
        return Tau(itertools.chain.from_iterable(blocks), self.ccomponent)


class ReducedTau:
    """ Tau in a reduction form """
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
    @property
    def Pzero(self) -> Iterable[Weight]:
        """ Search for weights w of V so that C_component + sum_k tau_red[w_k, k] = 0 """
        as_tau = Tau(self.values, self.ccomponent) # We use the scalar product of Tau
        for weight in Weight.all(self.small_d):
            if as_tau.dot_weight(weight) == 0:
                yield weight
