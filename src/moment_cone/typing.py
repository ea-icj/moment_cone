__all__ = (
    "Type", "NewType", "TypeVar", "TypeGuard", "Generic", "Optional",
    "TypeAlias", "TypeVarTuple", "Unpack",
    "cast", "overload",
    "Any", "Self", "Literal", "TYPE_CHECKING", "Concatenate",
    "ClassVar", "Final", "TypedDict",
    "Iterable", "Iterator", "Generator", "Mapping", "MutableMapping",
    "Sequence", "MutableSequence", "Sized", "Callable", "Hashable",
    "ABC", "abstractmethod",
    "T", "Ts", "T_co", "U",
    "Method",
)

from typing import Type, NewType, TypeVar, TypeGuard, Generic, Optional
from typing import TypeAlias, TypeVarTuple, Unpack
from typing import cast, overload
from typing import Any, Self, Literal, TYPE_CHECKING, Concatenate
from typing import ClassVar, Final, TypedDict
from collections.abc import Iterable, Iterator, Generator, Mapping, MutableMapping
from collections.abc import Sequence, MutableSequence, Sized, Callable, Hashable
from abc import ABC, abstractmethod

T = TypeVar('T') # Template parameter for generic classes or functions
#T_co = TypeVar('T_co', covariant=True)
Ts = TypeVarTuple("Ts")
T_co = TypeVar('T_co')
U = TypeVar('U') # Template parameter for generic classes or functions

# Method kind used in different parts of the algorithm
Method = Literal["probabilistic", "symbolic","symbolic_int"]

