from typing import Protocol, Generic, TypeAlias, ClassVar
from collections.abc import Iterable, Iterator

import numpy as np

from .typing import T

class SetOfVector(Protocol, Generic[T]):
    dtype: ClassVar[TypeAlias]
    def __init__(self, width: int): ...
    def add(self, element: Iterable[T]) -> bool: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> Iterator[Iterable[T]]: ...
    def clear(self) -> None: ...

class SetOfVector8(SetOfVector[np.int8]):
    dtype: ClassVar[TypeAlias] = np.int8
    ...

class SetOfVector16(SetOfVector[np.int16]):
    dtype: ClassVar[TypeAlias] = np.int16
    ...

class SetOfVector32(SetOfVector[np.int32]):
    dtype: ClassVar[TypeAlias] = np.int32
    ...

class SetOfVector64(SetOfVector[np.int64]):
    dtype: ClassVar[TypeAlias] = np.int64
    ...
