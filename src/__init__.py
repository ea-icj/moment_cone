from sage.all import SymmetricFunctions, QQ, Integer as SageInteger
sym_f = SymmetricFunctions(QQ).s()
import sage.libs.lrcalc.lrcalc as lr
import itertools
from dataclasses import dataclass

# Type hints
from typing import TypeVar, TypeAlias, Generic, Optional, Iterable, Iterator, cast
from collections.abc import Sequence
T = TypeVar('T')
Integer = TypeVar('Integer', int, SageInteger)
Weight: TypeAlias = Iterable[int]

def is_decreasing(l: Iterable[Integer]) -> bool:
    """ Check if a given sequence is not increasing """
    return all(a >= b for a, b in itertools.pairwise(l))

def compress(values: Iterable[T]) -> Iterable[tuple[T, int]]:
    """ Compress sequence of values by consecutive identical values and multiplicities """
    for value, group in itertools.groupby(values):
        yield value, sum(1 for _ in group)

def decompress(values: Iterable[T], mult: Iterable[int]) -> Iterable[T]:
    """ Decompress output from compress to the initial sequence """
    return itertools.chain.from_iterable(itertools.repeat(v, m) for v, m in zip(values, mult))

def trim_zeros(s: Sequence[Integer]) -> Sequence[Integer]:
    """ Remove trailing zeros from a sequence """
    for i, v in enumerate(reversed(s)):
        if v != 0:
            return s[:len(s) - i]
    else:
        return s
    
def all_weights(d: Iterable[int]) -> Iterable[Weight]:
    """ Returns all possible weights for a given sequence of length """
    return itertools.product(*(range(di) for di in d))


class PartialMatrix(Generic[T]):
    """ Kind of sparse matrix with fixed number of column, each of them containing a variable number of values """
    __slots__ = '_data', 'shape'
    _data: tuple[list[T], ...] # Internal representation of the columns
    shape: tuple[int, int] # Shape of the array (maximal number of rows and fixed number of cols)

    def __init__(self, nrows: int, ncols: int):
        """ Initialize empty partial matrix of given shape """
        self._data = tuple([] for _ in range(ncols))
        self.shape = nrows, ncols
    
    def size(self, idx: int) -> int:
        """ Number of rows for the given columns """
        return len(self._data[idx])
    
    @property
    def sizes(self) -> tuple[int, ...]:
        """ Number of rows for all columns """
        return tuple(len(col) for col in self._data)
    
    def __getitem__(self, idx: tuple[int, int]) -> T:
        """ Read-only access to an element using classical (i,j) order """
        i, j = idx
        return self._data[j][i]
    
    def column(self, j: int) -> Sequence[T]:
        """ Access to a given column """
        return self._data[j]
    
    def append(self, j: int, value: T) -> None:
        """ Append a value to the given column """
        self._data[j].append(value)

    def extend(self, j: int, values: Iterable[T]) -> None:
        """ Extend a column with given values """
        self._data[j].extend(values)

    @property
    def columns(self) -> Sequence[Sequence[T]]:
        """ Access to all columns """
        return self._data

    def __repr__(self) -> str:
        return " | ".join(" ".join(map(str, c)) for c in self.columns)


# TODO: un genre de FrozenPartialMatrix non mutable et mieux optimisé (on pourrait reprendre l'idée de tout mettre en ligne)

class Partition(Sequence[Integer]):
    """ Decreasing sequence of positive integer """
    __slots__ = "_data",
    _data: tuple[Integer, ...]

    def __init__(self, p: Sequence[Integer]):
        assert is_decreasing(p) and p[-1] >= 0, "Invalid partition"

        # Auto trim the partition
        self._data = tuple(trim_zeros(p))

    def __len__(self) -> int:
        return len(self._data)
    
    def __getitem__(self, idx: int) -> Integer:
        return self._data[idx] if idx < len(self) else cast(Integer, 0)
    
    def __iter__(self) -> Iterator[Integer]:
        return iter(self._data)
    
    def __repr__(self) -> str:
        return f"Partition({self._data})"


class Tau(Generic[Integer]):
    """ Tuple of partition along with a coefficient """
    __slots__ = 'ccomponent', '_components'
    ccomponent: Optional[Integer]
    _components: PartialMatrix[Integer]

    def __init__(self, components: Sequence[Sequence[Integer]] | PartialMatrix[Integer], ccomponent: Optional[Integer] = None):
        """ Tau initialization from a sequence of sub-group or directly from a partial matrix """
        self.ccomponent = ccomponent
        if isinstance(components, PartialMatrix):
            self._components = components # FIXME: copy ?
        else:
            self._components = PartialMatrix(max(len(c) for c in components), len(components))
            for j, c in enumerate(components):
                self._components.extend(j, c)
    
    def __len__(self) -> int:
        """ Number of components """
        return self._components.shape[1]

    @property
    def d(self) -> tuple[int, ...]:
        """ Length of each component """
        return self._components.sizes
    
    @property
    def components(self) -> Sequence[Sequence[Integer]]:
        """ Sequence of the components of tau """
        return self._components.columns

    @property
    def reduced(self) -> "ReducedTau[Integer]":
        """ Returns reduced form of tau """
        return ReducedTau(self)
    
    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + " | ".join(" ".join(map(str, c)) for c in self.components)
    
    def dot(self, weight: Weight) -> Integer:
        """ Scalar product of tau with a weight """
        assert self.ccomponent is not None
        return self.ccomponent, Integer + cast(Integer, sum(c[wi] for c, wi in zip(self.components, weight))) # type: ignore
    
    @property
    def is_regular(self) -> bool:
        """ Check is tau is regular assuming it is dominant """
        return not any(any(a == b for a, b in itertools.pairwise(c)) for c in self.components)
    
    def positive_weights(self, weights: Iterable[Weight]) -> dict[Integer, list[Weight]]:
        """ Inverse image of each non-zero p = <w, tau> for each w in weights """
        result = {}
        for chi in weights:
            p = self.dot(chi)
            if p >= 0:
                result.setdefault(p, []).append(chi)
        return result


class ReducedTau(Generic[Integer]):
    """ Tau in a reduction form """
    __slots__ = 'ccomponent', 'values', 'mult'
    ccomponent: Optional[Integer]
    values: PartialMatrix[Integer]
    mult: PartialMatrix[Integer]

    def __init__(self, tau: Tau[Integer]):

        values: PartialMatrix[Integer] = PartialMatrix(max(tau.d), len(tau))
        mult: PartialMatrix[Integer] = PartialMatrix(max(tau.d), len(tau))

        for j, component in enumerate(tau.components):
            for v, m in compress(component):
                values.append(j, v)
                mult.append(j, cast(Integer, m))

        self.ccomponent = tau.ccomponent
        self.values = values
        self.mult = mult

    def __len__(self) -> int:
        return self.values.shape[1]
    
    @property
    def small_d(self) -> tuple[int, ...]:
        return self.values.sizes
    
    def __getitem__(self, idx: tuple[int, int]) -> tuple[Integer, Integer]:
        return self.values[idx], self.mult[idx]
    
    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + " | ".join(
            " ".join(f"{v}^{m}" for v, m in zip(cv, cm))
            for cv, cm in zip(self.values.columns, self.mult.columns)
        )
    
    def Pzero(self, n: Integer) -> Iterable[Weight]:
        """ Search for weights w so that n + sum_k tau_red[w_k, k] = 0 """
        for weight in all_weights(self.small_d):
            s = cast(Integer, sum(self.values[i, j] for j, i in enumerate(weight)))
            if s + n == 0:
                yield weight

    

def fct_Weights_of_Nu(nu: PartialMatrix[Partition[Integer]]) -> PartialMatrix[Integer]:
    """ Constructing size_nu from nu """
    size_nu = PartialMatrix[Integer](*nu.shape)

    for j, col in enumerate(nu.columns):
        size_nu.extend(j, (cast(Integer, sum(p)) for p in col))

    return size_nu



##### Pour 444, il y a un cas où on rentre dans LR_multi avec L=[] (liste vide) et nu=[].
## Cela bugue mais je ne comprends pas si c'est normal d'avoir une telle entrée
