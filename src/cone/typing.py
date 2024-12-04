from typing import NewType, TypeVar, Generic, Optional, Iterable, Iterator, TypeGuard
from typing import cast, overload
from typing import Any, Self, TYPE_CHECKING
from collections.abc import Sequence, MutableSequence, Sized, Callable

T = TypeVar('T') # Template parameter for generic classes or functions
