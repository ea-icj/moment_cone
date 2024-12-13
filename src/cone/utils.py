from .typing import *
import itertools
import functools
import operator
import copy as cp
from sympy.utilities.iterables import multiset_permutations

__all__ = (
    "is_decreasing",
    "is_increasing",
    "group_by_block",
    "expand_blocks",
    "trim_zeros",
    "count",
    "prod",
    "short_prod",
    "Embeddings_mod_sym",
    "extend_with_repetitions",
    "flatten_dictionary",
    "grading_dictionary",
    "Are_Isom_Mod",
    "Is_Sub_Mod",
    "quotient_C_Mod",
    "dictionary_list_lengths",
    
)

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
    """ Classical product of all given values """
    return functools.reduce(operator.mul, values)

def short_prod(values: Iterable[int]) -> int:
    """ Product of value with sort-circuit if result is 0 """
    result = 0
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

def filter_dict_by_key(d: dict[T, U], predicate: Callable[[T], bool]) -> dict[T, U]:
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

def Embeddings_mod_sym(d: Sequence[int], e: Sequence[int])-> Iterable["Permutation"]:
    """
    List of permutations of e that are at most d

    d and e are list of integers of the same length  (typically, dimensions), each in decreasing order.

    Returns the list of permutation of e (each encoded by a permutation of the indices) such that the value in the i-th component of the permuted e is at most d[i]
    Outputs are irredundant modulo symmetries of e and d
    
    Example:
    >>> d = [4, 4, 3, 3, 2]
    >>> e = [4, 3, 3, 2, 1]
    >>> emb = list(Embeddings_mod_sym(d, e))
    >>> for pe in emb:
    ...     print(pe)
    Permutation((0, 1, 2, 3, 4))
    Permutation((0, 1, 2, 4, 3))
    Permutation((0, 3, 1, 2, 4))
    Permutation((0, 4, 1, 2, 3))
    >>> for pe in emb:
    ...     print(pe(e))
    (4, 3, 3, 2, 1)
    (4, 3, 3, 1, 2)
    (4, 2, 3, 3, 1)
    (4, 1, 3, 3, 2)
    """
    # TODO: lowercase name, move to Permutation
    from .permutation import Permutation

    eg = list(group_by_block(e))
    dg = list(group_by_block(d))
    partial_sum_mult_d = [0] + list(itertools.accumulate([x for _,x in dg]))
    partial_sum_mult_e = [0] + list(itertools.accumulate([x for _,x in eg]))
    indices_eg = expand_blocks([i for i, x in enumerate(eg)], [x[1] for i, x in enumerate(eg)]) # same as e, but with repetition of indices, rather than values
    for ep in multiset_permutations(indices_eg):
        p_i = cp.copy(partial_sum_mult_e)[:-1] # FIXME: slice of list ix already a copy
        indices_e = []
        for i in range(len(e)):
           indices_e.append(p_i[ep[i]])
           p_i[ep[i]] += 1
        if all(e[indices_e[i]] <= d[i] for i in range(len(e))) and all(is_increasing(indices_e[a:b]) for a,b in itertools.pairwise(partial_sum_mult_d)):
           yield Permutation(indices_e)

def extend_with_repetitions(seq:Sequence[T], l: int) -> Iterable[tuple[T, ...]]:
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
    
def flatten_dictionary(dic: dict[U, Iterable[T]]) -> list[T]:
    """
    Returns the concatenation of all list stored as values in a dict.
    
    Example:
    >>> d = {0: [1, 2], 1: [4, 5], 2: [3, 6, 7]}
    >>> sorted(flatten_dictionary(d))
    [1, 2, 3, 4, 5, 6, 7]
    """
    # TODO: can we return the iterable directly?
    return list(itertools.chain.from_iterable(dic.values()))
   
def dictionary_list_lengths(dic: dict[U, Sequence[T]]) -> dict[U, int]:
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
   
def Is_Sub_Mod(M1: dict[int, int], M2: dict[int, int]) -> bool:
    """
    Kind of order on dictionary of int -> int.

    Examples:
    >>> d1 = {0: 0, 1: 1, 2: 2}
    >>> d2 = {0: 0, 2: 2}
    >>> Is_Sub_Mod(d1, d2)
    False
    >>> d3 = {0: 0, 1: 1, 2: 1}
    >>> Is_Sub_Mod(d1, d3)
    False
    >>> d4 = {0: 0, 1: 1, 2: 3}
    >>> Is_Sub_Mod(d1, d4)
    True
    """
    # TODO: lowercase name, move to more specific file?
    for p in M1.keys():
        if p not in M2.keys() or M1[p] > M2[p] :
            return False
    return True

def Are_Isom_Mod(M1 : dict[int, int], M2 : dict[int, int]) -> bool:
    """ Comparison of dictionary of int -> int.
    
    Examples:
    >>> d1 = {0: 0, 1: 1, 2: 2}
    >>> d2 = {0: 0, 2: 2}
    >>> Are_Isom_Mod(d1, d2)
    False
    >>> d3 = {0: 0, 1: 1, 2: 1}
    >>> Are_Isom_Mod(d1, d3)
    False
    >>> d4 = {0: 0, 1: 1, 2: 3}
    >>> Are_Isom_Mod(d1, d4)
    False
    >>> d5 = {0: 0, 1: 1, 2: 2}
    >>> Are_Isom_Mod(d1, d5)
    True
    """
    # TODO: lowercase name, move to more specific file, probably useless since it is simply a comparison
    return M1 == M2
    
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


  

