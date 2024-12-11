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
    "extend_with_repetitions"
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
    


def Embeddings_mod_sym(d: Sequence[int],e:Sequence[int])-> Iterable[Iterable[int]]: # List of permutations of e that are at most d
    """
    d and e are list of integers the same length  (typically, dimensions), each in decreasing order 
    returns the list of permutation of e (each encoded by a permutation of the indices) such that the in the i-th component of the permuted e is at most d[i]
    Outputs are irredundant modulo symmetries of e and d
    
    Example:
    >>>d=[4,4,3,3,2]
    >>>e=[4,3,3,2,1]
    >>>emb=Embeddings(d,e)
    >>>emb
    [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 3, 1, 2, 4], [0, 4, 1, 2, 3]]
    >>>[e[emb[-1][i]] for i in range(5)]
    [4, 1, 3, 3, 2]
    """
    Res=[]
    eg=list(group_by_block(e))
    dg=list(group_by_block(d))
    partial_sum_mult_d=[0]+list(itertools.accumulate([x for _,x in dg]))
    partial_sum_mult_e=[0]+list(itertools.accumulate([x for _,x in eg]))
    indices_eg=expand_blocks([i for i,x in enumerate(eg)],[x[1] for i,x in enumerate(eg)]) #same as e, but with repetition of indices, rather than values
    for ep in multiset_permutations(indices_eg):
        p_i=cp.copy(partial_sum_mult_e)[:-1]
        indices_e=[]
        for i in range(len(e)):
           indices_e.append(p_i[ep[i]])
           p_i[ep[i]]+=1
        #Res.append(indices_e)
    #return(Res)
        if all(e[indices_e[i]] <= d[i] for i in range(len(e))) and all(is_increasing(indices_e[a:b]) for a,b in itertools.pairwise(partial_sum_mult_d)):
           Res.append(indices_e) # Une chance sur deux inverser ep et e si nécessaire 
    return(Res)

def extend_with_repetitions(seq:Sequence[T],l:[int])->Sequence[T]:
    """from a sequence seq of length <=l with no repetition, returns the list of all expanded sequences of length l obtained from seq by repetitions of some elements."""
    if len(seq)==1 : 
       return([l*seq])
    if len(seq)==l : 
       return([seq])
    Res=[]
    for i in range(l-len(seq)+1):
        Res+=[(i+1)*[seq[0]]+ext for ext in extend_with_repetitions(seq[1:],l-i-1)]
    return(Res)

