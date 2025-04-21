import contextlib
from multiprocessing.pool import Pool
import itertools

from .typing import *

class Parallel(contextlib.AbstractContextManager["Parallel"]):
    """
    Parallel computation context with some convenient algorithms (map, filter, ...)

    If the number of jobs is 1, then these algorithms are essentially calls to built-in functions.

    Note that passes functions must be pickle and thus declared in the top level of a module
    (not a local function or a lambda).

    Example:
    >>> with Parallel(1) as p:
    ...     r = p.map(round, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ...     r = list(r) # r in only an iterable
    >>> r
    [0, 0, 0, 0, 0, 1, 1, 1, 1]

    >>> from operator import gt # >= operator
    >>> with Parallel(6) as p:
    ...     r = p.filter(gt, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], 0.4)
    ...     r = list(r) # results must be retrieved before the end of the parallel context
    >>> r
    [0.5, 0.6, 0.7, 0.8, 0.9]
    """
    pool: Optional[Pool]

    def __init__(self, n_jobs: Optional[int] = None):
        """ Starting a parallel computation context
        
        If n_jobs is None, it lets multiprocessing.pool.Pool choose the number of processes
        depending on your CPU and the execution context.
        """
        assert n_jobs is None or n_jobs > 0

        if n_jobs is None or n_jobs > 1:
            self.pool = Pool(n_jobs)
        else:
            self.pool = None

    def __enter__(self) -> Self:
        """ Entering context """
        return self
    
    def terminate(self) -> None:
        """ Terminating all processes """
        if self.pool is not None:
            self.pool.close()
            self.pool.terminate()

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        """ Leaving context """
        self.terminate()

    @staticmethod
    def _starmap(args: tuple[Callable[[Unpack[Ts]], T], Unpack[Ts]]) -> T:
        return args[0](*args[1:])

    def map(self, func: Callable[Concatenate[Any, ...], T], /, *iterables: Iterable[Any], chunk_size: int = 1, unordered: bool = False) -> Iterable[T]:
        """
        Like the built-in map, it applies given function with one parameter from each iterable

        By default, results are returned in the same order as the input
        but they can be returned when ready when setting unordered to True.
        """
        if self.pool is None:
            return map(func, *iterables)
        else:
            if unordered:
                imap = self.pool.imap_unordered
            else:
                imap = self.pool.imap

            return imap(
                Parallel._starmap,
                zip(itertools.repeat(func), *iterables),
                chunksize=chunk_size,
            )

    @staticmethod
    def _filter(args: tuple[Callable[[T, Unpack[Ts]], U], T, Unpack[Ts]]) -> tuple[U, T]:
        return args[0](*args[1:]), args[1]
    
    def filter(self, func: Callable[Concatenate[Any, ...], bool], iterable: Iterable[T], /, *args: Any, chunk_size: int = 1, unordered: bool = False) -> Generator[T]:
        """
        Like the built-in filter, it returns the elements from iterable for whose func returns True.

        The element is passed as the first parameter to func, the additional args (not iterables)
        are passed as supplementary parameters to each call of func.

        By default, results are returned in the same order as the input
        but they can be returned when ready when setting unordered to True.
        """
        from itertools import repeat
        if self.pool is None:
            for el in iterable:
                if func(el, *args):
                    yield el
        else:
            if unordered:
                imap = self.pool.imap_unordered
            else:
                imap = self.pool.imap

            results = imap(
                Parallel._filter,
                zip(repeat(func), iterable, *(repeat(a) for a in args)),
                chunksize=chunk_size,
            )
            for keep, el in results:
                if keep:
                    yield el
    