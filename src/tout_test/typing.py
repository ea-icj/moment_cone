__all__ = (
    "Type", "NewType", "TypeVar", "TypeGuard", "Generic", "Optional", "TypeAlias",
    "cast", "overload",
    "Any", "Self", "Literal", "TYPE_CHECKING", "Concatenate",
    "ClassVar",
    "Iterable", "Iterator", "Generator", "Mapping",
    "Sequence", "MutableSequence", "Sized", "Callable",
    "T", "T_co", "U",
    "Method",
)

from typing import Type, NewType, TypeVar, TypeGuard, Generic, Optional, TypeAlias
from typing import cast, overload
from typing import Any, Self, Literal, TYPE_CHECKING, Concatenate
from typing import ClassVar
from collections.abc import Iterable, Iterator, Generator, Mapping
from collections.abc import Sequence, MutableSequence, Sized, Callable

T = TypeVar('T') # Template parameter for generic classes or functions
#T_co = TypeVar('T_co', covariant=True)
T_co = TypeVar('T_co')
U = TypeVar('U') # Template parameter for generic classes or functions

# Method kind used in different parts of the algorithm
Method = Literal["probabilistic", "symbolic"]

