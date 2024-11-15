from .typing import *
import itertools

def is_decreasing(l: Iterable[int]) -> bool:
    """ Check if a given sequence is not increasing """
    return all(a >= b for a, b in itertools.pairwise(l))

def compress(values: Iterable[T]) -> Iterable[tuple[T, int]]:
    """ Compress sequence of values by consecutive identical values and multiplicities """
    for value, group in itertools.groupby(values):
        yield value, sum(1 for _ in group)

def decompress(values: Iterable[T], mult: Iterable[int]) -> Iterable[T]:
    """ Decompress output from compress to the initial sequence """
    return itertools.chain.from_iterable(itertools.repeat(v, m) for v, m in zip(values, mult))

def trim_zeros(s: Sequence[int]) -> Sequence[int]:
    """ Remove trailing zeros from a sequence """
    for i, v in enumerate(reversed(s)):
        if v != 0:
            return s[:len(s) - i]
    else:
        return ()
    
