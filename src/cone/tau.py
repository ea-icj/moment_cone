from .typing import *
from .matrix import PartialMatrix
from .weight import Weight
from .root import Root

import itertools
from functools import cached_property

"""
TODO:
- test file
- same method for both scalar product using @singledispatch?
- property or method for Pzero?
"""

class Tau:
    """ Tuple of partition along with a coefficient """
    #__slots__ = 'ccomponent', '_components' # FIXME cached_property cannot work without __dict__ ... => self managed cache or removing __slots__
    ccomponent: Optional[int]
    _components: PartialMatrix[int]

    def __init__(self, components: Sequence[Sequence[int]] | PartialMatrix[int], ccomponent: Optional[int] = None):
        """ Tau initialization from a sequence of sub-group or directly from a partial matrix """
        self.ccomponent = ccomponent
        if isinstance(components, PartialMatrix):
            self._components = components # FIXME: copy ?
        else:
            self._components = PartialMatrix(max(len(c) for c in components), len(components))
            for j, c in enumerate(components):
                self._components.extend(j, c)

    @staticmethod
    def from_flatten(s: Iterable[int], d: Dimension) -> "Tau":
        """ Returns tau from its flattened version """
        all_components = tuple(s)
        sum_dim = sum(d)
        if len(all_components) == sum_dim:
            ccomponent = None
        elif len(all_components) == sum_dim + 1:
            ccomponent = all_components[0]
            all_components = all_components[1:]
        else:
            raise ValueError("Invalid number of components")
        
        # TODO: from_flatten in PartialMatrix so that it can be optimized
        # accordingly to the internal storage strategy!
        shift = itertools.accumulate(d, initial=0)
        return Tau(
            tuple(all_components[a:b] for a, b in itertools.pairwise(shift)),
            ccomponent
        )
    
    def __len__(self) -> int:
        """ Number of components """
        return self._components.shape[1]

    @cached_property
    def d(self) -> Dimension:
        """ Length of each component """
        return self._components.sizes
    
    @cached_property
    def components(self) -> tuple[list[int], ...]:
        """ Sequence of the components of tau """
        return self._components.columns
    
    @property
    def flattened(self) -> Iterable[int]:
        """
        Returns the whole tau as a unique sequence.
        
        ie (cc-component | column1 | column2 | ...)
        """
        if self.ccomponent is None:
            return itertools.chain.from_iterable(self.components)
        else:
            return itertools.chain((self.ccomponent,), *self.components)

    @cached_property
    def reduced(self) -> "ReducedTau":
        """ Returns reduced form of tau """
        return ReducedTau(self)
    
    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + " | ".join(" ".join(map(str, c)) for c in self.components)
    
    # TODO: or same dot method with @singledispatch ?
    def dot_weight(self, weight: Weight, ccomponent: Optional[int] = None) -> int:
        """ Scalar product of tau with a weight of V """
        assert not (self.ccomponent is None and ccomponent is None)
        if ccomponent is None:
            ccomponent = self.ccomponent
        return cast(int, ccomponent) + sum(c[wi] for c, wi in zip(self.components, weight))
    
    def dot_root(self, root: Root) -> int:
        """ Scalar product of tau with a weight of U """
        c = self.components[root.k]
        return c[root.i] - c[root.j]
    
    @cached_property
    def is_regular(self) -> bool:
        """ Check if tau is regular assuming it is dominant """
        # FIXME: why not just all(a > b for ...) then?
        return not any(any(a == b for a, b in itertools.pairwise(c)) for c in self.components)
    
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

    def positive_weights(self, weights: Iterable[Weight]) -> dict[int, list[Weight]]:
        """ Inverse image of each non-zero p = <w, tau> for each w in weights """
        result = {}
        for chi in weights:
            p = self.dot_weight(chi)
            if p >= 0:
                result.setdefault(p, []).append(chi)
        return result
    
    @cached_property
    def sort_mod_sym_dim(self) -> "Tau":
        """ Sort tau by block of the dimensions """
        from .utils import compress
        col_split = itertools.accumulate((length for _, length in compress(self.d)), initial=0)
        columns = itertools.chain.from_iterable(
            sorted(self.components[a:b])
            for a, b in itertools.pairwise(col_split)
        )
        return Tau(tuple(columns), self.ccomponent)


class ReducedTau:
    """ Tau in a reduction form """
    __slots__ = 'ccomponent', 'values', 'mult'
    ccomponent: Optional[int]
    values: PartialMatrix[int]
    mult: PartialMatrix[int]

    def __init__(self, tau: Tau):
        from .utils import compress
        values: PartialMatrix[int] = PartialMatrix(max(tau.d), len(tau))
        mult: PartialMatrix[int] = PartialMatrix(max(tau.d), len(tau))

        for j, component in enumerate(tau.components):
            for v, m in compress(component):
                values.append(j, v)
                mult.append(j, cast(int, m))

        self.ccomponent = tau.ccomponent
        self.values = values
        self.mult = mult

    def __len__(self) -> int:
        return self.values.shape[1]
    
    @cached_property
    def small_d(self) -> Dimension:
        return self.values.sizes
    
    def __getitem__(self, idx: tuple[int, int]) -> tuple[int, int]:
        return self.values[idx], self.mult[idx]
    
    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + " | ".join(
            " ".join(f"{v}^{m}" for v, m in zip(cv, cm))
            for cv, cm in zip(self.values.columns, self.mult.columns)
        )
    
    @property # FIXME: returns an generator and not a sequence. Is it adapted to a property (or a cached_property?)
    def Pzero(self) -> Iterable[Weight]:
        """ Search for weights w of V so that C_component + sum_k tau_red[w_k, k] = 0 """
        as_tau = Tau(self.values, self.ccomponent) # We use the scalar product of Tau
        for weight in Weight.all(self.small_d):
            if as_tau.dot_weight(weight) == 0:
                yield weight