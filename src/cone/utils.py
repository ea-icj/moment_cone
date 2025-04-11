from .typing import *

import itertools
import functools
import operator
import copy as cp

__all__ = (
    "is_decreasing",
    "is_increasing",
    "group_by_block",
    "expand_blocks",
    "trim_zeros",
    "count",
    "prod",
    "short_prod",
    "extend_with_repetitions",
    "flatten_dictionary",
    "grading_dictionary",
    "compare_C_Mod",
    "quotient_C_Mod",
    "dictionary_list_lengths",
    "unique_combinations",
    "to_literal",
    "get_function_by_name",
    "line_profiler",
    "CachedClass",
    "fl_dic",
)

if TYPE_CHECKING:
    from line_profiler import LineProfiler
    import pstats


def is_decreasing(l: Iterable[int]) -> bool:
    """ Check if a given sequence is not increasing """
    return all(a >= b for a, b in itertools.pairwise(l))

def is_increasing(l: Iterable[int]) -> bool:
    """ Check if a given sequence is not decreasing """
    return all(a <= b for a, b in itertools.pairwise(l))

def group_by_block(values: Iterable[T]) -> Iterable[tuple[T, int]]:
    """
    Compress sequence of values by consecutive identical values and multiplicities
    
    For each block, returns a pair of (value, multiplicity).
    """
    for value, group in itertools.groupby(values):
        yield value, sum(1 for _ in group)

def expand_blocks(values: Iterable[T], mult: Iterable[int]) -> Iterable[T]:
    """ Decompress output from `group_by_block` to the initial sequence """
    return itertools.chain.from_iterable(itertools.repeat(v, m) for v, m in zip(values, mult))

def symmetries(values: Iterable[T]) -> Iterable[int]:
    """
    Returns the multiplicities of identical consecutive values in a sequence
    
    Example:
    >>> tuple(symmetries((4, 4, 3, 2, 2, 2, 1)))
    (2, 1, 3, 1)
    """
    return (length for _, length in group_by_block(values))

def trim_zeros(s: Sequence[int]) -> Sequence[int]:
    """ Remove trailing zeros from a sequence """
    for i, v in enumerate(reversed(s)):
        if v != 0:
            return s[:len(s) - i]
    else:
        return ()
    
def count(s: Iterable[T]) -> int:
    """ Count number of elements in an iterable """
    if isinstance(s, Sized):
        return len(s)
    else:
        return sum(1 for _ in s) # It seems that they exist faster method using `collections.deque`

def prod(values: Iterable[int]) -> int:
    """
    Classical product of all given values
    
    >>> prod([1, 2, 3])
    6
    >>> prod([1, 2, 3, 0, 5])
    0
    """
    return functools.reduce(operator.mul, values)

def short_prod(values: Iterable[int]) -> int:
    """
    Product of value with sort-circuit if result is 0
    
    >>> short_prod([1, 2, 3])
    6
    >>> short_prod([1, 2, 3, 0, 5])
    0

    Testing short-circuit feature:
    >>> g = iter([1, 2, 3, 0, 5])
    >>> short_prod(g)
    0
    >>> next(g)
    5
    """
    result = 1
    for v in values:
        result *= v
        if result == 0:
            return 0
    return result


def grading_dictionary(elements: Iterable[T], fn: Callable[[T], U]) -> dict[U, list[T]]:
    """ 
    From a sequence of elements and a function that applies on these elements,
    generates a dictionary that maps each image to it's preimage.

    Example:
    >>> elements = range(40)
    >>> fn = lambda e: e % 7
    >>> gd = grading_dictionary(elements, fn)
    >>> for k in sorted(gd.keys()):
    ...     print(f"{k}:", gd[k])
    0: [0, 7, 14, 21, 28, 35]
    1: [1, 8, 15, 22, 29, 36]
    2: [2, 9, 16, 23, 30, 37]
    3: [3, 10, 17, 24, 31, 38]
    4: [4, 11, 18, 25, 32, 39]
    5: [5, 12, 19, 26, 33]
    6: [6, 13, 20, 27, 34]
    """
    result: dict[U, list[T]] = {}
    for e in elements:
        v = fn(e)
        result.setdefault(v, []).append(e)
    return result

def filter_dict_by_key(d: Mapping[T, U], predicate: Callable[[T], bool]) -> dict[T, U]:
    """
    Filter a dictionary using a predicate on its keys
    Example:
    >>> elements = range(40)
    >>> fn = lambda e: e % 7
    >>> gd = grading_dictionary(elements, fn)
    >>> gdf = filter_dict_by_key(gd, lambda k: k % 2 == 0)
    >>> for k in sorted(gdf.keys()):
    ...     print(f"{k}:", gdf[k])
    0: [0, 7, 14, 21, 28, 35]
    2: [2, 9, 16, 23, 30, 37]
    4: [4, 11, 18, 25, 32, 39]
    6: [6, 13, 20, 27, 34]
    """
    return {k: v for k, v in d.items() if predicate(k)}

def extend_with_repetitions(seq: Sequence[T], l: int) -> Iterable[tuple[T, ...]]:
    """
    From a sequence seq of length <= l with no repetition, returns the list of all expanded sequences of length l obtained from seq by repetitions of some elements.

    Examples:
    >>> for l in extend_with_repetitions([1, 2, 3], 5):
    ...     print(l)
    (1, 2, 3, 3, 3)
    (1, 2, 2, 3, 3)
    (1, 2, 2, 2, 3)
    (1, 1, 2, 3, 3)
    (1, 1, 2, 2, 3)
    (1, 1, 1, 2, 3)
    >>> for l in extend_with_repetitions([1], 5):
    ...     print(l)
    (1, 1, 1, 1, 1)
    """
    assert 0 < len(seq) <= l, "Incompatible sequence length and/or target length"
    if len(seq) == 1: 
       yield l * (seq[0],)
    elif len(seq) == l: 
       yield tuple(seq)
    else:
        for i in range(l - len(seq) + 1):
            for tail in extend_with_repetitions(seq[1:], l - i - 1):
                yield (i + 1) * (seq[0],) + tail
    
def flatten_dictionary(dic: Mapping[U, Iterable[T]]) -> list[T]:
    """
    Returns the concatenation of all list stored as values in a dict.
    
    Example:
    >>> d = {0: [1, 2], 1: [4, 5], 2: [3, 6, 7]}
    >>> sorted(flatten_dictionary(d))
    [1, 2, 3, 4, 5, 6, 7]
    """
    # TODO: can we return the iterable directly?
    return list(itertools.chain.from_iterable(dic.values()))
   
def dictionary_list_lengths(dic: Mapping[U, Sequence[T]]) -> dict[U, int]:
    """
    From a dictionary of list, returns the dictionary of the length of each list.

    Example:
    >>> d = {0: [1, 2], 1: [4, 5], 2: [3, 6, 7]}
    >>> dl = dictionary_list_lengths(d)
    >>> for k in sorted(dl.keys()):
    ...     print((k, dl[k]))
    (0, 2)
    (1, 2)
    (2, 3)
    """
    return {key: len(value) for key, value in dic.items()}
   

def compare_C_Mod(
        M1: dict[int, int],
        M2: dict[int, int],
        relation: Callable[[int, int], bool]
    ) -> bool :
    # TODO: lowercase, docstring, move to list_of_W ?
    return all(
        relation(M1.get(key, 0), M2.get(key, 0))
        for key in set(M1) | set(M2)
    )
    

def quotient_C_Mod(M1 : dict[int, int], M2 : dict[int, int]) -> dict[int, int]:
    """ Quotient of two dictionary int -> int.

    Examples:
    >>> d1 = {0: 0, 1: 1, 2: 2}
    >>> d2 = {0: 0, 2: 2}
    >>> quotient_C_Mod(d1, d2)
    {1: 1}
    >>> d3 = {0: 0, 1: 1, 2: 1}
    >>> quotient_C_Mod(d1, d3)
    {2: 1}
    >>> d4 = {0: 1, 1: 1, 2: 3}
    >>> quotient_C_Mod(d1, d4)
    {0: -1, 2: -1}
    >>> d5 = {0: 0, 1: 1, 2: 2}
    >>> quotient_C_Mod(d1, d5)
    {}
    """
    # TODO: lowercase name, move to more specific file?
    M: dict[int, int] = {}
    for p in M1.keys():
        # Default value for M2[p] : 0
        difference = M1[p] - M2.get(p, 0)
        # Add only non-zero differences
        if difference != 0:
            M[p] = difference
    return M

def multiset_permutations(m: Iterable[T]) -> Generator[list[T]]:
    """
    Returns the unique permutations of the given multiset m
    
    Is it simply a wrapper of sympy.utilities.iterables.multiset_permutations
    so that to get the correct return type.
    """
    from sympy.utilities.iterables import multiset_permutations as mp
    return cast(Generator[list[T]], mp(m)) # type: ignore

def orbit_symmetries(flatten: Iterable[T], symmetries: Iterable[int]) -> Generator[Iterable[T]]:
    """
    Permutation inside each block of given sizes

    Note that each returned permutation is only an iterable (not a list or tuple).
    
    If this is too slow, we may consider the remarks/propositions from:
    - https://stackoverflow.com/questions/19676109/how-to-generate-all-the-permutations-of-a-multiset/
    - https://stackoverflow.com/questions/70057504/speed-up-multiset-permutations

    Example:
    >>> orbits = orbit_symmetries((2, 2, 4, 1, 2, 1, 4), (3, 3, 1))
    >>> for p in orbits:
    ...     print(tuple(p))
    (2, 2, 4, 1, 1, 2, 4)
    (2, 2, 4, 1, 2, 1, 4)
    (2, 2, 4, 2, 1, 1, 4)
    (2, 4, 2, 1, 1, 2, 4)
    (2, 4, 2, 1, 2, 1, 4)
    (2, 4, 2, 2, 1, 1, 4)
    (4, 2, 2, 1, 1, 2, 4)
    (4, 2, 2, 1, 2, 1, 4)
    (4, 2, 2, 2, 1, 1, 4)
    """
    from .blocks import Blocks
    blocks = (multiset_permutations(block) for block in Blocks.from_flatten(tuple(flatten), symmetries))
    for p in itertools.product(*blocks):
        yield itertools.chain.from_iterable(p)
  

def unique_combinations(mylist: Sequence[int], k: int) -> list[tuple[int, ...]]:
    """
    The list is sorted.
    The list is viewed as a multiset. 
    Create the list of multisets contained in list of cardinality k

    Example :
    >>> unique_combinations([3,3,3,2,2,1],3)
    [(3, 3, 3), (3, 3, 2), (3, 3, 1), (3, 2, 2), (3, 2, 1), (2, 2, 1)]
    """
    def backtrack(start: int, current: list[int]) -> None:
        if len(current) == k:
            result.append(tuple(current))
            return
        for i in range(start, len(mylist)):
            # Éviter de répéter un élément au même niveau
            if i > start and mylist[i] == mylist[i - 1]:
                continue
            current.append(mylist[i])
            backtrack(i + 1, current)
            current.pop()

    result: list[tuple[int, ...]] = []
    backtrack(0, [])
    return result


def to_literal(literals: Any, arg: str) -> str:
    """ Convert a given string to the corresponding literal
    
    Returns first literal that begins with given string.
    """
    from typing import get_args, get_origin
    assert get_origin(literals) == Literal
    for s in get_args(literals):
        assert isinstance(s, str)
        if s.lower().startswith(arg.lower()):
            return s

    raise ValueError(f"No literal matches {arg}")


def get_function_by_name(name: str) -> Callable[..., Any]:
    """ Returns a function by it's name 
    
    The name may include full module or a sub-module of cone
    and can be a method of class.
    """
    # Trying in global scope
    try:
        return cast(Callable[..., Any], globals()[name])
    except KeyError:
        pass

    # Trying in objects imported from cone
    import cone
    try:
        return cast(Callable[..., Any], getattr(cone, name))
    except AttributeError:
        pass

    # Trying with sub-module or class name
    from importlib import import_module
    for full_name in (name, "cone." + name):
        components = full_name.split(".")
        # Splitting module name and object name
        for n in range(1, len(components)):
            # Loading first part as a module
            try:
                mod = import_module(".".join(components[:n]))
            except ModuleNotFoundError:
                break

            # Get remaining part as multi-level objects (eg class)
            try:
                obj = getattr(mod, components[n])
                for cname in components[n+1:]:
                    obj = getattr(obj, cname)
                return cast(Callable[..., Any], obj)
            except AttributeError:
                pass
    
    raise ValueError(f"No function {name} found")


def line_profiler(
        profiled_fn: Iterable[str | Callable[..., Any]] | str | Callable[..., Any],
        called_fn: Callable[..., T],
        *args: Any,
        **kwargs: Any,
        ) -> tuple[T, "LineProfiler"]:
    """
    Profile given functions by line when launching called_fn with given arguments
    
    It returns the result of the called_fn and the LineProfiler instance.
    Statistics can be printed by calling the print_stats method of the returned
    LineProfiler.
    """
    if isinstance(profiled_fn, str) or callable(profiled_fn):
        profiled_fn = [profiled_fn]

    from line_profiler import LineProfiler
    lp = LineProfiler()
    for name_or_fn in profiled_fn:
        lp.add_function(
            get_function_by_name(name_or_fn)
            if isinstance(name_or_fn, str)
            else name_or_fn
        )

    result: T = lp.runcall(called_fn, *args, **kwargs) # type: ignore
    return result, lp


def cprofile(
        called_fn: Callable[..., T],
        *args: Any,
        file_name: Optional[str] = None,
        **kwargs: Any,
        ) -> tuple[T, "pstats.Stats"]:
    """
    Profile calls and hierarchy of calls of the functions executed during
    computation of given called_fn.
    
    Beside the result of the the called_fn, it also returns a instance
    of pstats.Stats and can also output to Stats and kcachegrind format
    if a file_name is given.
    """
    from cProfile import Profile
    from pstats import Stats
    with Profile() as prof:
        result = called_fn(*args, **kwargs)
        stats = Stats(prof)

    if file_name is not None:
        stats.dump_stats(file_name + ".pstats")
        from pyprof2calltree import convert # type: ignore
        convert(stats, file_name + ".prof")

    return result, stats

def fl_dic(D : dict[int,list["Root"]],L: list[int])-> dict[int,int]:
    """
    dictionnary h -> sum lenght D[i] for i <=h, for h in L.
    L has to be a sorted list.
    """
    length_D = {}
    total = 0
    for key in sorted(D.keys()) :
        total += len(D[key])
        length_D[key] = total
    result = {}
    for j in L:
        # Chercher la plus grande clé i ≤ j (ou 0 si aucune clé ne convient)
        i_max = max((i for i in sorted(D.keys()) if i <= j), default=None)
        result[j] = length_D[i_max] if i_max is not None else 0
    return(result)    
    

class CachedClass:
    """ Base class that ensure instance uniqueness relatively to the construction arguments 
    
    It relies on a private class dictionary that maps the construction arguments
    to a corresponding instance. In order to work, all arguments must be hashable.

    Remark that when mixing positional and keyword arguments, it may leads to
    different instances of a same object depending on the order and the type
    of the arguments.
    """
    __all_instances: ClassVar[
        dict[tuple[tuple[Any, ...], tuple[tuple[str, Any], ...]], Self]
    ] = {}

    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        return cls.__all_instances.setdefault(
            (args, tuple(kwargs.items())),
            super().__new__(cls),
        )

