"""
Overlay on Numpy array in 1D and 2D so that to get better type annotation,
especially for array with optional value.
"""

__all__ = (
    "Array1D", "Array2D", "Array3D",
    "EmptyArray1D", "EmptyArray2D", "EmptyArray3D", "empty_array",
    "OptionalArray1D", "OptionalArray2D", "OptionalArray3D", "optional_array",
)

import numpy as np
from numpy.typing import NDArray, ArrayLike

from .typing import *


class Array1D(Generic[T_co]):
    data: NDArray[Any]
    
    def __init__(self, data: Iterable[T_co], dtype: Optional[Type[T_co]] = None, **kwargs: Any):
        self.data = np.array(data, dtype=dtype, **kwargs)
        assert self.data.ndim == 1

    @property
    def shape(self) -> tuple[int]:
        return cast(tuple[int], self.data.shape)

    @property
    def dtype(self) -> Type[T_co]:
        return cast(Type[T_co], self.data.dtype)
    
    def copy(self) -> "Array1D[T_co]":
        return Array1D(self.data.copy(), dtype=self.dtype)

    @overload
    def __getitem__(self, i: int) -> T_co:
        ...

    @overload
    def __getitem__(self, i: slice) -> "Array1D[T_co]":
        ...

    def __getitem__(self, i: int | slice) -> T_co | "Array1D[T_co]":
        if isinstance(i, int):
            return cast(T_co, self.data[i])
        else:
            return Array1D(self.data[i])

    def __iter__(self) -> Iterator[T_co]:
        return iter(self.data)

    def __len__(self) -> int:
        return len(self.data)
    
    def tolist(self) -> list[T_co]:
        return cast(list[T_co], self.data.tolist())
            
    def __repr__(self) -> str:
        return f"Array1D({self.data})"


class Array2D(Generic[T_co]):
    data: NDArray[Any]
    
    def __init__(self, data: Iterable[Iterable[T_co]], dtype: Optional[Type[T_co]] = None, **kwargs: Any):
        self.data = np.array(data, dtype=dtype, **kwargs)
        assert self.data.ndim == 2

    @property
    def shape(self) -> tuple[int, int]:
        return cast(tuple[int, int], self.data.shape)

    @property
    def dtype(self) -> Type[T_co]:
        return cast(Type[T_co], self.data.dtype)
    
    def copy(self) -> "Array2D[T_co]":
        return Array2D(self.data.copy(), dtype=self.dtype)
    
    def ravel(self) -> Array1D[T_co]:
        return Array1D(self.data.ravel())
    
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
            return cast(T_co, self.data[i, j])
        elif isinstance(i, slice) and isinstance(j, slice):
            return Array2D(self.data[i, j])
        else:
            return Array1D(self.data[i, j])
        
    def __iter__(self) -> Iterator[Array1D[T_co]]:
        return iter(self.data)
    
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
        self.data[*idx] = data

    def __repr__(self) -> str:
        return f"Array2D({self.data})"
    

class Array3D(Generic[T_co]):
    data: NDArray[Any]
    
    def __init__(self, data: Iterable[Iterable[T_co]], dtype: Optional[Type[T_co]] = None, **kwargs: Any):
        self.data = np.array(data, dtype=dtype, **kwargs)
        assert self.data.ndim == 3

    @property
    def shape(self) -> tuple[int, int, int]:
        return cast(tuple[int, int, int], self.data.shape)

    @property
    def dtype(self) -> Type[T_co]:
        return cast(Type[T_co], self.data.dtype)
    
    def copy(self) -> "Array3D[T_co]":
        return Array3D(self.data.copy(), dtype=self.dtype)
        
    def ravel(self) -> Array1D[T_co]:
        return Array1D(self.data.ravel())
    
    @overload
    def __getitem__(self, idx: tuple[int, int, int]) -> T_co:
        ...

    @overload
    def __getitem__(self, idx: tuple[int, int, slice] | tuple[int, slice, int] | tuple[slice, int, int]) -> Array1D[T_co]:
        ...

    @overload
    def __getitem__(self, idx: int | tuple[int, slice, slice] | tuple[slice, int, slice] | tuple[slice, slice, int]) -> "Array2D[T_co]":
        ...

    @overload
    def __getitem__(self, idx: slice | tuple[slice, slice, slice]) -> "Array3D[T_co]":
        ...

    def __getitem__(self, idx: int | slice | tuple[int | slice, int | slice, int | slice]) -> T_co | Array1D[T_co] | "Array2D[T_co]" | "Array3D[T_co]":
        if isinstance(idx, tuple):
            i, j, k = idx
        else:
            i, j, k = idx, slice(None, None, None), slice(None, None, None)

        int_cnt = sum(1 for _ in (i, j, k) if isinstance(_, int))
        match int_cnt:
            case 0: return Array3D(self.data[i, j, k])
            case 1: return Array2D(self.data[i, j, k])
            case 2: return Array1D(self.data[i, j, k])
            case 3: return cast(T_co, self.data[i, j, k])
            case _: raise ValueError()
        
    def __iter__(self) -> Iterator[Array1D[T_co]]:
        return iter(self.data)
    
    @overload
    def __setitem__(self, idx: tuple[int, int, int], data: T_co) -> None:
        ...

    @overload
    def __setitem__(self, idx: tuple[int, int, slice] | tuple[int, slice, int] | tuple[slice, int, int], data: Iterable[T_co]) -> None:
        ...
    
    @overload
    def __setitem__(self, idx: tuple[int, slice, slice] | tuple[slice, int, slice] | tuple[slice, slice, int], data: Iterable[Iterable[T_co]]) -> None:
        ...

    @overload
    def __setitem__(self, idx: tuple[slice, slice, slice], data: Iterable[Iterable[Iterable[T_co]]]) -> None:
        ...

    def __setitem__(self, idx: tuple[int | slice, int | slice, int | slice ], data: T_co | Iterable[T_co] | Iterable[Iterable[T_co]] | Iterable[Iterable[Iterable[T_co]]]) -> None:
        self.data[*idx] = data

    def __repr__(self) -> str:
        return f"Array2D({self.data})"
    

OptionalArray1D: TypeAlias = Array1D[Optional[T_co]]
OptionalArray2D: TypeAlias = Array2D[Optional[T_co]]
OptionalArray3D: TypeAlias = Array3D[Optional[T_co]]

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

@overload
def optional_array(
        shape: tuple[int, int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> OptionalArray3D[T_co]:
    ...

def optional_array(
        shape: int | tuple[int] | tuple[int, int] | tuple[int, int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> OptionalArray1D[T_co] | OptionalArray2D[T_co] | OptionalArray3D[T_co]:
    """ Creates an array (1D or 2D only) with optional value """
    data = np.full(shape, None, dtype=object, **kwargs)
    if isinstance(shape, tuple) and len(shape) > 1:
        match len(shape):
            case 2: return Array2D[Optional[T_co]](data)
            case 3: return Array3D[Optional[T_co]](data)
            case _: raise NotImplementedError()
    else:
        return Array1D[Optional[T_co]](data)



EmptyArray1D: TypeAlias = Array1D[T_co]
EmptyArray2D: TypeAlias = Array2D[T_co]
EmptyArray3D: TypeAlias = Array3D[T_co]

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

@overload
def empty_array(
        shape: tuple[int, int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> EmptyArray3D[T_co]:
    ...

def empty_array(
        shape: int | tuple[int] | tuple[int, int] | tuple[int, int, int],
        dtype: Type[T_co],
        **kwargs: Any
        ) -> EmptyArray1D[T_co] | EmptyArray2D[T_co] | EmptyArray3D[T_co]:
    """ Creates an uninitialized array (1D or 2D only) """
    data = np.empty(shape, dtype=dtype, **kwargs)
    if isinstance(shape, tuple) and len(shape) > 1:
        match len(shape):
            case 2: return Array2D[T_co](data)
            case 3: return Array3D[T_co](data)
            case _: raise NotImplementedError()
    else:
        return Array1D[T_co](data)




