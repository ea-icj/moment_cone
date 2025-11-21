from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

from .linear_group import LinearGroup
from .tau import Tau
from .typing import *


class UniqueTau(ABC):
    """ Archetype of unique storage of tau """
    G: LinearGroup

    def __init__(self, G: LinearGroup):
        self.G = G
    
    @abstractmethod
    def add(self, tau: Tau) -> bool:
        """ Add a tau to the set
        
        Returns True if the wasn't already in the set.
        """
        ...

    def __call__(self, tau: Tau) -> bool:
        return self.add(tau)
    
    @abstractmethod
    def __iter__(self) -> Iterator[Tau]:
        """ Iterate over all unique tau """
        ...

    @abstractmethod
    def __len__(self) -> int:
        """ Number of unique tau """
        ...

    @abstractmethod
    def clear(self) -> None:
        """ Clear stored tau (to help reducing memory usage) """
        ...

    def __repr__(self) -> str:
        return f"{type(self).__name__}(#tau={len(self)})"

    @staticmethod
    def _find_best_dtype(tau: Tau, dtypes: Iterable[type[np.integer]]) -> tuple[int, type[np.integer]]:
        """ Helper method to find best suitable integer representation """
        min_coeff = min(tau.flattened)
        max_coeff = max(tau.flattened)
        for i, dtype in enumerate(dtypes):
            dtype_info = np.iinfo(dtype)
            if dtype_info.min <= min_coeff and max_coeff <= dtype_info.max:
                return i, dtype
        raise ValueError("Unsupported integer width")


class SetOfTau(UniqueTau):
    """ Storage of unique tau based on a simple Python set of Tau
    

    Example:

    >>> G = LinearGroup((3, 2))
    >>> sot = SetOfTau(G)
    >>> print(sot)
    SetOfTau(#tau=0)
    >>> sot.add(Tau.from_flatten((1, 2, 3, 4, 1032), G))
    True
    >>> sot.add(Tau.from_flatten((1, 2, 2, 4, -2456), G))
    True
    >>> sot.add(Tau.from_flatten((1, 2, 3, 4, 1032), G))
    False
    >>> print(sot)
    SetOfTau(#tau=2)
    >>> for tau in sot:
    ...     print(tau)
    1 2 3 | 4 1032
    1 2 2 | 4 -2456
    >>> sot.clear()
    >>> print(sot)
    SetOfTau(#tau=0)
    """
    data: set[Tau]

    def __init__(self, G: LinearGroup) -> None:
        super().__init__(G)
        self.data = set()

    def __len__(self) -> int:
        return len(self.data)
    
    def __iter__(self) -> Iterator[Tau]:
        return iter(self.data)
    
    def add(self, tau: Tau) -> bool:
        if tau not in self.data:
            self.data.add(tau)
            return True
        else:
            return False
        
    def clear(self) -> None:
        self.data.clear()


DType = TypeVar("DType",bound=np.generic)
class NumpyVector(Generic[DType]):
    """
    Class that mimic a C++ vector for tabular data using Numpy
    
    Example :

    >>> import numpy as np
    >>> v = NumpyVector(3, np.int8)
    >>> v.size
    0
    >>> v.push_back([1, 2, 3])
    >>> v.size
    1
    >>> v[0]
    array([1, 2, 3], dtype=int8)
    >>> v.push_back([2, 3, 4])
    >>> v[1]
    array([2, 3, 4], dtype=int8)
    >>> len(v)
    2
    >>> v.clear()
    >>> len(v)
    0
    """
    shape: tuple[int, ...] #! shape of each element
    _data: NDArray[DType]
    dtype: type[DType]
    size: int

    def __init__(self, shape: int | tuple[int, ...], dtype: type[DType], *args: Any, **kwargs: Any):
        import numpy as np
        if isinstance(shape, int):
            shape = (shape,)
        self.shape = shape
        self._data = np.empty((0,) + self.shape, dtype, *args, **kwargs)
        self.dtype = dtype
        self.size = 0

    @property
    def capacity(self) -> int:
        return self._data.shape[0]

    def reserve(self, size: int) -> None:
        assert size >= self.size
        self._data.resize((size,) + self.shape, refcheck=False)

    def __getitem__(self, idx: int) -> NDArray[DType]:
        assert idx < self.size
        return self._data[idx, ...]
    
    def __len__(self) -> int:
        return self.size
    
    def push_back(self, element: Iterable[DType] | NDArray[DType]) -> None:
        element = np.asarray(element)
        assert element.shape == self.shape
        if self.size == self.capacity:
            self.reserve(max(1, self.capacity * 2))
        self._data[self.size, ...] = element
        self.size += 1
    
    def pop_back(self) -> None:
        assert self.size > 0
        self.size -=1

    def insert(self, element: Iterable[DType] | NDArray[DType], idx: int) -> None:
        element = np.asarray(element)
        assert element.shape == self.shape
        assert idx <= self.size
        if self.size == self.capacity:
            self.reserve(max(1, self.capacity * 2))
        self._data[idx+1:self.size+1] = self._data[idx:self.size]
        self._data[idx] = element
        self.size += 1
        
    def __iter__(self) -> Iterator[NDArray[DType]]:
        return iter(self._data[:self.size])
    
    def __repr__(self) -> str:
        return f"NumpyVector(shape={self.shape}, dtype={self.dtype.__name__}, data={self._data[:self.size]})"
    
    def clear(self) -> None:
        self._data.resize((0,) + self.shape, refcheck=False)
        self.size = 0


class MinimalSetOfTau(UniqueTau):
    """ Memory optimized set of Tau by using a set of minimal elements
    
    and optimized integer representation.

    Example:

    >>> G = LinearGroup((3, 2))
    >>> sot = MinimalSetOfTau(G)
    >>> print(sot)
    MinimalSetOfTau(#tau=0)
    >>> sot.add(Tau.from_flatten((1, 2, 3, 4, 1032), G))
    True
    >>> sot.add(Tau.from_flatten((1, 2, 2, 4, -2456), G))
    True
    >>> sot.add(Tau.from_flatten((1, 2, 3, 4, 1032), G))
    False
    >>> print(sot)
    MinimalSetOfTau(#tau=2)
    >>> for tau in sot:
    ...     print(tau)
    1 2 3 | 4 1032
    1 2 2 | 4 -2456
    >>> sot.clear()
    >>> print(sot)
    MinimalSetOfTau(#tau=0)
    """
    _data: tuple[NumpyVector[np.generic], ...]
    _sets: tuple[set[Any], ...]
    _element_type: list[type]

    def __init__(
            self,
            G: LinearGroup,
            dtypes: tuple[type[np.integer], ...] = (np.int8, np.int16, np.int32, np.int64),
        ):
        super().__init__(G)
        self.G = G
        self._data = tuple(
            NumpyVector(sum(G), dtype)
            for dtype in dtypes + (np.object_,)
        )
        self._sets = tuple(set() for dtype in dtypes + (np.object_,))
        self._element_type = []

        # Custom set element based on index in each Numpy vector
        from dataclasses import dataclass
        for data in self._data:
            @dataclass(slots=True)
            class DataElement:
                data: ClassVar[NumpyVector[Any]]
                idx: int
                def __hash__(self) -> int:
                    return hash(tuple(self.data[int(self.idx)]))
                
                def __eq__(self, other: object) -> bool:
                    if not isinstance(other, type(self)):
                        return NotImplemented
                    return cast(bool, np.all(self.data[int(self.idx)] == other.data[int(other.idx)]))

            DataElement.data = data
            self._element_type.append(DataElement)
    
    def add(self, tau: Tau) -> bool:
        """ Add a tau to the set
        
        Returns True if the wasn't already in the set.
        """
        # The idea is to add tau as a new element in the storage
        # then to check if it was already there (using a DataElement as key)
        # and to remove it from the storage if it is not new
        try:
            best_dtype, _ = self._find_best_dtype(
                tau,
                (cast(type[np.integer], v.dtype) for v in self._data[:-1])
            )
        except ValueError:
            best_dtype = len(self._data) - 1 # Python int
        data = self._data[best_dtype]

        data.push_back(np.asarray(tau.flattened))
    
        element = self._element_type[best_dtype](len(data) - 1)
        found = element in self._sets[best_dtype]
        if found:
            data.pop_back()
        else:
            self._sets[best_dtype].add(element)
        return not found

    def __iter__(self) -> Iterator[Tau]:
        from itertools import chain
        for el in chain.from_iterable(self._data):
            yield Tau.from_flatten(map(int, el), self.G)

    def __len__(self) -> int:
        return sum(len(d) for d in self._data)
    
    def clear(self) -> None:
        for d, s in zip(self._data, self._sets):
            d.clear()
            s.clear()


class SetOfTauCpp(UniqueTau):
    """ Memory optimized set of Tau by using a C++ std::set and minimal elements
     
    and also optimized integer representations.

    Example:

    >>> G = LinearGroup((3, 2))
    >>> sot = SetOfTauCpp(G)
    >>> print(sot)
    SetOfTauCpp(#tau=0)
    >>> sot.add(Tau.from_flatten((1, 2, 3, 4, 1032), G))
    True
    >>> sot.add(Tau.from_flatten((1, 2, 2, 4, -2456), G))
    True
    >>> sot.add(Tau.from_flatten((1, 2, 3, 4, 1032), G))
    False
    >>> print(sot)
    SetOfTauCpp(#tau=2)
    >>> for tau in sot:
    ...     print(tau)
    1 2 3 | 4 1032
    1 2 2 | 4 -2456
    >>> sot.clear()
    >>> print(sot)
    SetOfTauCpp(#tau=0)
    """
    from . import _set_of_vector as sov
    G: LinearGroup
    _data: tuple[sov.SetOfVector8, sov.SetOfVector16, sov.SetOfVector32, sov.SetOfVector64] 

    def __init__(
            self,
            G: LinearGroup,
        ):
        from . import _set_of_vector as sov
        self.G = G
        width = sum(G)
        self._data = sov.SetOfVector8(width), sov.SetOfVector16(width), sov.SetOfVector32(width), sov.SetOfVector64(width)
    
    def add(self, tau: Tau) -> bool:
        """ Add a tau to the set
        
        Returns True if the wasn't already in the set.
        If drill is True, then no tau is added even if not already in the set.
        """
        # The idea is to add tau as a new element in the storage
        # then to check if it was already there (using an Element as key)
        # and to remove it from the storage if it is not new
        dtype_idx, _ = self._find_best_dtype(tau, (np.int8, np.int16, np.int32, np.int64))
        data = self._data[dtype_idx]
        added = data.add(tau.flattened) # type: ignore
        return added

    def clear(self) -> None:
        for data in self._data:
            data.clear()

    def __iter__(self) -> Iterator[Tau]:
        from itertools import chain
        for el in chain.from_iterable(self._data):
            yield Tau.from_flatten(map(int, el), self.G) #type: ignore

    def __len__(self) -> int:
        return sum(len(d) for d in self._data)
    

UniqueTauStr = Literal[
    "SetOfTau",
    "MinimalSetOfTau",
    "SetOfTauCpp",
]

unique_tau_dict: Final[dict[UniqueTauStr, type[UniqueTau]]] = {
    "SetOfTau": SetOfTau,
    "MinimalSetOfTau": MinimalSetOfTau,
    "SetOfTauCpp": SetOfTauCpp,
}
