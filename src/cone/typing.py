from typing import NewType, TypeVar, Generic, Optional, Iterable, Iterator, TypeGuard
from typing import cast, overload
from typing import Any, Self, Literal, TYPE_CHECKING
from collections.abc import Sequence, MutableSequence, Sized, Callable

T = TypeVar('T') # Template parameter for generic classes or functions

U = TypeVar('U') # Template parameter for generic classes or functions

# Method kind used in different parts of the algorithm
Method = Literal["probabilistic", "symbolic"]

