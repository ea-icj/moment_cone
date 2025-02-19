"""
Overlay on Numpy array in 1D and 2D so that to get better type annotation,
especially for array with optional value.
"""

__all__ = (
    "Array1D", "Array2D",
    "EmptyArray1D", "EmptyArray2D", "empty_array",
    "OptionalArray1D", "OptionalArray2D", "optional_array",
)

import numpy as np
from numpy.typing import NDArray, ArrayLike

from .typing import *


class Array1D(Generic[T_co]):
    __data: NDArray[Any]
    
    def __init__(self, data: Iterable[T_co], dtype: Optional[Type[T_co]] = None, **kwargs: Any):
        self.__data = np.array(data, dtype=dtype, **kwargs)
        assert self.__data.ndim == 1

    @property
    def shape(self) -> tuple[int]:
        return cast(tuple[int], self.__data.shape)
    
    @overload
    def __getitem__(self, i: int) -> T_co:
        ...

    @overload
    def __getitem__(self, i: slice) -> "Array1D[T_co]":
        ...

    def __getitem__(self, i: int | slice) -> T_co | "Array1D[T_co]":
        if isinstance(i, int):
            return cast(T_co, self.__data[i])
        else:
            return Array1D(self.__data[i])

    def __iter__(self) -> Iterator[T_co]:
        return iter(self.__data)

    def __len__(self) -> int:
        return len(self.__data)
    
    def tolist(self) -> list[T_co]:
        return cast(list[T_co], self.__data.tolist())
    
    def __repr__(self) -> str:
        return f"Array1D({self.__data})"


class Array2D(Generic[T_co]):
    __data: NDArray[Any]
    
    def __init__(self, data: Iterable[Iterable[T_co]], dtype: Optional[Type[T_co]] = None, **kwargs: Any):
        self.__data = np.array(data, dtype=dtype, **kwargs)
        assert self.__data.ndim == 2

    @property
    def shape(self) -> tuple[int, int]:
        return cast(tuple[int, int], self.__data.shape)
    
    def ravel(self) -> Array1D[T_co]:
        return Array1D(self.__data.ravel())
    
    @overload
    def __getitem__(self, idx: tuple[int, int]) -> T_co:
        ...

    @overload
    def __getitem__(self, idx: int | tuple[int, slice]) -> Array1D[T_co]:
        ...

    @overload
    def __getitem__(self, idx: tuple[slice, int]) -> Array1D[T_co]:
        ...

    @overload
    def __getitem__(self, idx: tuple[slice, slice]) -> "Array2D[T_co]":
        ...

    def __getitem__(self, idx: int | slice | tuple[int | slice, int | slice]) -> T_co | Array1D[T_co] | "Array2D[T_co]":
        if isinstance(idx, tuple):
            i, j = idx
        else:
            i, j = idx, slice(None, None, None)

        if isinstance(i, int) and isinstance(j, int):
            return cast(T_co, self.__data[i, j])
        elif isinstance(i, slice) and isinstance(j, slice):
            return Array2D(self.__data[i, j])
        else:
            return Array1D(self.__data[i, j])
        
    def __iter__(self) -> Iterator[Array1D[T_co]]:
        return iter(self.__data)
    
    @overload
    def __setitem__(self, idx: tuple[int, int], data: T_co) -> None:
        ...

    @overload
    def __setitem__(self, idx: tuple[int, slice] | tuple[slice, int], data: Iterable[T_co]) -> None:
        ...
    
    @overload
    def __setitem__(self, idx: tuple[slice, slice], data: Iterable[Iterable[T_co]]) -> None:
        ...

    def __setitem__(self, idx: tuple[int | slice, int | slice], data: T_co | Iterable[T_co] | Iterable[Iterable[T_co]]) -> None:
        self.__data[*idx] = data

    def __repr__(self) -> str:
        return f"Array2D({self.__data})"
    

OptionalArray1D: TypeAlias = Array1D[Optional[T_co]]
OptionalArray2D: TypeAlias = Array2D[Optional[T_co]]

@overload
def optional_array(
        shape: int | tuple[int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> OptionalArray1D[T_co]:
    ...

@overload
def optional_array(
        shape: tuple[int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> OptionalArray2D[T_co]:
    ...

def optional_array(
        shape: int | tuple[int] | tuple[int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> OptionalArray1D[T_co] | OptionalArray2D[T_co]:
    """ Creates an array (1D or 2D only) with optional value """
    data = np.full(shape, None, dtype=object, **kwargs)
    if isinstance(shape, tuple) and len(shape) > 1:
        return Array2D[Optional[T_co]](data)
    else:
        return Array1D[Optional[T_co]](data)



EmptyArray1D: TypeAlias = Array1D[T_co]
EmptyArray2D: TypeAlias = Array2D[T_co]

@overload
def empty_array(
        shape: int | tuple[int],
        dtype: Type[T_co],
        **kwargs: Any,
        ) -> EmptyArray1D[T_co]:
    ...

@overload
def empty_array(
        shape: tuple[int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> EmptyArray2D[T_co]:
    ...

def empty_array(
        shape: int | tuple[int] | tuple[int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> EmptyArray1D[T_co] | EmptyArray2D[T_co]:
    """ Creates an uninitialized array (1D or 2D only) """
    data = np.empty(shape, dtype=object, **kwargs)
    if isinstance(shape, tuple) and len(shape) > 1:
        return Array2D[T_co](data)
    else:
        return Array1D[T_co](data)




