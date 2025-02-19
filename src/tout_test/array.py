"""
Overlay on Numpy array in 1D and 2D so that to get better type annotation
"""

__all__ = (
    "Array1D", "Array2D",
    "EmptyArray1D", "EmptyArray2D",
)

import numpy as np
from numpy.typing import NDArray, ArrayLike

from .typing import *

class Array1D(Generic[T]):
    __data: NDArray
    
    def __init__(self, data: Iterable[T], *args, dtype: Optional[Type[T]] = None, **kwargs):
        self.__data = np.array(data, *args, dtype=dtype, **kwargs)

    @property
    def shape(self) -> tuple[int]:
        return cast(tuple[int], self.__data.shape)
    
    @overload
    def __getitem__(self, i: int) -> T:
        ...

    @overload
    def __getitem__(self, i: slice) -> "Array1D[T]":
        ...

    def __getitem__(self, i: int | slice) -> T | "Array1D[T]":
        if isinstance(i, int):
            return self.__data[i]
        else:
            return Array1D(self.__data[i])

    def __iter__(self) -> Iterator[T]:
        return iter(self.__data)

    def __len__(self) -> int:
        return len(self.__data)
    
    def tolist(self) -> list[T]:
        return self.__data.tolist()


class Array2D(Generic[T]):
    __data: NDArray
    
    def __init__(self, data: Iterable[Iterable[T]], *args, dtype: Optional[Type[T]] = None, **kwargs):
        self.__data = np.array(data, *args, dtype=dtype, **kwargs)

    @property
    def shape(self) -> tuple[int, int]:
        return cast(tuple[int, int], self.__data.shape)
    
    def ravel(self) -> Array1D[T]:
        return Array1D(self.__data.ravel())
    
    @overload
    def __getitem__(self, idx: tuple[int, int]) -> T:
        ...

    @overload
    def __getitem__(self, idx: int | tuple[int, slice]) -> Array1D[T]:
        ...

    @overload
    def __getitem__(self, idx: tuple[slice, int]) -> Array1D[T]:
        ...

    @overload
    def __getitem__(self, idx: tuple[slice, slice]) -> "Array2D[T]":
        ...

    def __getitem__(self, idx: int | slice | tuple[int | slice, int | slice]) -> T | Array1D[T] | "Array2D[T]":
        if isinstance(idx, tuple):
            i, j = idx
        else:
            i, j = idx, slice(None, None, None)

        if isinstance(i, int) and isinstance(j, int):
            return self.__data[i, j]
        elif isinstance(i, slice) and isinstance(j, slice):
            return Array2D(self.__data[i, j])
        else:
            return Array1D(self.__data[i, j])
        
    def __iter__(self) -> Iterator[Array1D]:
        return iter(self.__data)
    
    @overload
    def __setitem__(self, idx: tuple[int, int], data: T) -> None:
        ...

    @overload
    def __setitem__(self, idx: tuple[int, slice] | tuple[slice, int], data: Iterable[T]) -> None:
        ...
    
    @overload
    def __setitem__(self, idx: tuple[slice, slice], data: Iterable[Iterable[T]]) -> None:
        ...

    def __setitem__(self, idx: tuple[int | slice, int | slice], data: T | Iterable[T] | Iterable[Iterable[T]]) -> None:
        self.__data[*idx] = data


def EmptyArray1D(shape: int, dtype: Type[T], **kwargs) -> Array1D[T]:
    return Array1D[T](np.empty(shape, dtype=dtype, **kwargs))

def EmptyArray2D(shape: tuple[int, int], dtype: Type[T], **kwargs) -> Array2D[T]:
    return Array2D[T](np.empty(shape, dtype=dtype, **kwargs))

def test1(data: Iterable[T]):
    print(data)

def test2(data: Sequence[T]):
    print(data)

def test3(data: Iterable[Iterable[T]]):
    print(data)

def toto():
    a = Array1D([], dtype=int)
    b = EmptyArray1D(3, dtype=int)
    aa = a[0]
    bb = a[:1]
    c = Array2D([[1, 2, 3], [4, 5, 6]], dtype=int)
    cc = c[0]
    dd = c[:, 1]
    c[0, :] = [1, 2]
    test3(c)

    for i, v in zip(a, c):
        w = cast(Array1D[int], v)
        
