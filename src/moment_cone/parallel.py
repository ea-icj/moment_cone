__all__ = (
    "Parallel",
    "ParallelExecutorStr",
)

from contextlib import AbstractContextManager
import concurrent.futures as futures
from multiprocessing.pool import Pool
from multiprocessing import Manager, Queue
from multiprocessing.managers import SyncManager
from abc import ABC
import itertools
from argparse import ArgumentParser, Namespace

from .typing import *



class ParallelExecutor(AbstractContextManager["ParallelExecutor"], ABC):
    """
    Base class for a parallel (or sequential) task executor
    
    Support manual creation/termination or through a context.
    """
    max_workers: Optional[int]
    chunk_size: int
    unordered: bool

    def __init__(
            self, max_workers: Optional[int] = None,
            chunk_size: int = 1,
            unordered: bool = False,
        ):
        """ Default parameters for an executor
        
        max_workers: maximal number of tasks that run in parallel
        chunk_size: maximal number of tasks in the queue of each worker

        These parameters can be ignored for some executor (eg Sequential)
        """
        self.max_workers = max_workers
        self.chunk_size = chunk_size
        self.unordered = unordered

    @property
    def is_parallel(self) -> bool:
        return True

    def shutdown(self, wait: bool = True) -> None:
        """ Shutdown executor and possibly wait that all tasks are finished """
        pass

    def __enter__(self) -> Self:
        """ Entering context """
        return self
    
    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        """ Leaving context """
        self.shutdown()
    
    @abstractmethod
    def map(self, 
            fn: Callable[..., T],
            /,
            *iterables: Iterable[Any],
            chunk_size: Optional[int] = None,
            unordered: Optional[bool] = None) -> Iterable[T]:
        """ Launch fn on each parameters from the given iterables
        
        This method should yield the results when available (like multiprocessing.Pool.imap),
        possibly out of order if unordered is True.

        If chunk_size is not given, the one given at the executor creation is used.
        """
        ...

    @staticmethod
    def _filter(args: tuple[Callable[[T, Unpack[Ts]], U], T, Unpack[Ts]]) -> tuple[U, T]:
        """ Internal method to ease filtering """
        return args[0](*args[1:]), args[1]
    
    def filter(self,
               fn: Callable[[T, Unpack[Ts]], bool],
               iterable: Iterable[T],
               /,
               *args: Unpack[Ts],
               chunk_size: Optional[int] = None,
               unordered: Optional[bool] = None) -> Iterable[T]:
        """ Filter the given iterable by fn
        
        Additional arguments may be given through args. Results may be yielded
        out of order if unordered is True.

        If chunk_size is not given, the one given at the executor creation is used.
        """
        from itertools import repeat
        results = self.map(
            ParallelExecutor._filter,
            zip(repeat(fn), iterable, *(repeat(a) for a in args)),
            chunk_size=chunk_size,
            unordered=unordered
        )
        for keep, value in results:
            if keep:
                yield value


class SequentialExecutor(ParallelExecutor):
    """ Sequential executor (rely on map) """  
    @property
    def is_parallel(self) -> bool:
        return False

    def map(self, 
            fn: Callable[Concatenate[Any, ...], T],
            /,
            *iterables: Iterable[Any],
            chunk_size: Optional[int] = None,
            unordered: Optional[bool] = None) -> Iterable[T]:
        return map(fn, *iterables)
    

class MultiProcessingPoolExecutor(ParallelExecutor):
    """
    Parallel computation context with some convenient algorithms (map, filter, ...)

    If the number of jobs is 1, then these algorithms are essentially calls to built-in functions.

    Note that passes functions must be pickle and thus declared in the top level of a module
    (not a local function or a lambda).

    Example:
    
    >>> with MultiProcessingPoolExecutor(1) as p:
    ...     r = p.map(round, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ...     r = list(r) # r in only an iterable
    >>> r
    [0, 0, 0, 0, 0, 1, 1, 1, 1]

    >>> from operator import gt # >= operator
    >>> with MultiProcessingPoolExecutor(6) as p:
    ...     r = p.filter(gt, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], 0.4)
    ...     r = list(r) # results must be retrieved before the end of the parallel context
    >>> r
    [0.5, 0.6, 0.7, 0.8, 0.9]
    """
    pool: "Pool"

    def __init__(self,
                max_workers: Optional[int] = None,
                chunk_size: int = 1,
                unordered: bool = False,
                *args: Any,
                **kwargs: Any):
        """ Starting a parallel computation context
        
        If n_jobs is None, it lets multiprocessing.pool.Pool choose the number of processes
        depending on your CPU and the execution context.
        """
        from multiprocessing import Pool
        super().__init__(
            max_workers=max_workers,
            chunk_size=chunk_size,
            unordered=unordered
        )
        self.pool = Pool(self.max_workers, *args, **kwargs)
    
    def shutdown(self, wait: bool = True) -> None:
        """ Terminating all processes """
        pool = self.pool
        pool.close()
        if wait:
            pool.join()
        else:
            pool.terminate()

    @staticmethod
    def _starmap(args: tuple[Callable[[Unpack[Ts]], T], Unpack[Ts]]) -> T:
        return args[0](*args[1:])

    def map(self,
            func: Callable[Concatenate[Any, ...], T],
            /, 
            *iterables: Iterable[Any],
            chunk_size: Optional[int] = None,
            unordered: Optional[bool] = None) -> Iterable[T]:
        """
        Like the built-in map, it applies given function with one parameter from each iterable

        By default, results are returned in the same order as the input
        but they can be returned when ready when setting unordered to True.
        """
        chunk_size = chunk_size or self.chunk_size
        if unordered is None: unordered = self.unordered
        if unordered:
            imap = self.pool.imap_unordered
        else:
            imap = self.pool.imap

        return imap(
            MultiProcessingPoolExecutor._starmap,
            zip(itertools.repeat(func), *iterables),
            chunksize=chunk_size,
        )


class FutureParallelExecutor(ParallelExecutor, ABC):
    """ Parallel executor based on the concurrent.futures module """
    executor: futures.Executor

    def __init__(self, executor: futures.Executor, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self.executor = executor

    def shutdown(self, wait: bool = True) -> None:
        self.executor.shutdown(wait=wait, cancel_futures=not wait)

    @staticmethod
    def chunk_task(fn: Callable[Concatenate[Any, ...], T],
                   chunk_args: Iterable[tuple[Any, ...]],
                   ) -> list[T]:
        return [fn(*args) for args in chunk_args]

    def map(self, 
            fn: Callable[Concatenate[Any, ...], T],
            /,
            *iterables: Iterable[Any],
            chunk_size: Optional[int] = None,
            unordered: Optional[bool] = None) -> Iterable[T]:
        chunk_size = chunk_size or self.chunk_size
        if unordered is None: unordered = self.unordered
        if not unordered:
            yield from self.executor.map(fn, *iterables, chunksize=chunk_size)
        else:
            """
            Previous implementation was :
                all_futures = (self.executor.submit(fn, *args) for args in zip(*iterables))
            with the use of futures.as_completed but it seems to force the complete
            evaluation of all_futures before continuing, thus breaking the lazy features.

            The implementation below relies on futures.wait instead and replace finished
            tasks by new ones to keep a pool of running tasks.

            As before, chunk_size is not used and performance may be worse but the lazy
            computation now works.
            """
            from concurrent.futures import Future
            from itertools import islice
            not_done: set[Future[list[T]]] = set()
            all_args = iter(zip(*iterables))
            max_workers: int = self.executor._max_workers # type: ignore
            for _ in range(max_workers):
                chunk_args = list(islice(all_args, chunk_size))
                if not chunk_args: break
                not_done.add(self.executor.submit(
                    FutureParallelExecutor.chunk_task,
                    fn,
                    chunk_args,
                ))

            while not_done:
                done, not_done = futures.wait(not_done, return_when=futures.FIRST_COMPLETED)
                for future in done:
                    yield from future.result()
                    chunk_args = list(islice(all_args, chunk_size))
                    if not chunk_args: continue
                    not_done.add(self.executor.submit(
                        FutureParallelExecutor.chunk_task,
                        fn,
                        chunk_args,
                    ))


class FutureProcessExecutor(FutureParallelExecutor):
    def __init__(self,
                max_workers: Optional[int] = None,
                chunk_size: int = 1,
                unordered: bool = False,
                *args: Any,
                **kwargs: Any):
        executor = futures.ProcessPoolExecutor(max_workers, *args, **kwargs)
        super().__init__(
            executor,
            max_workers=max_workers,
            chunk_size=chunk_size,
            unordered=unordered,
        )


class FutureThreadExecutor(FutureParallelExecutor):
    def __init__(self,
                max_workers: Optional[int] = None,
                chunk_size: int = 1,
                unordered: bool = False,
                *args: Any,
                **kwargs: Any):
        executor = futures.ThreadPoolExecutor(max_workers, *args, **kwargs)
        super().__init__(
            executor,
            max_workers=max_workers,
            chunk_size=chunk_size,
            unordered=unordered,
        )


class FutureMPIExecutor(FutureParallelExecutor):
    def __init__(self,
                max_workers: Optional[int] = None,
                chunk_size: int = 1,
                unordered: bool = False,
                *args: Any,
                **kwargs: Any):
        from mpi4py.futures import MPIPoolExecutor
        executor = MPIPoolExecutor(max_workers, *args, **kwargs)
        super().__init__(
            executor,
            max_workers=max_workers,
            chunk_size=chunk_size,
            unordered=unordered,
        )


class MultiProcessingQueueExecutor(ParallelExecutor):
    """ Trying to parallelize using input and output queues
    
    Not very efficient in practice...
    """
    manager: SyncManager

    def __init__(
            self,
            max_workers: Optional[int] = None,
            chunk_size: int = 1,
            unordered: bool = False,
        ):
        if max_workers is None:
            from multiprocessing import cpu_count
            max_workers = cpu_count()
        super().__init__(max_workers, chunk_size, unordered)
        self.manager = Manager()

    def shutdown(self, wait: bool = True):
        self.manager.shutdown()

    def map(self, 
            fn: Callable[..., T],
            /,
            *iterables: Iterable[Any],
            chunk_size: Optional[int] = None,
            unordered: Optional[bool] = None) -> Iterable[T]:
        chunk_size = chunk_size or self.chunk_size
        if unordered is None: unordered = self.unordered

        # Initializing queues
        assert self.max_workers is not None
        max_workers = self.max_workers
        print("max_workers =", max_workers)
        queue_max_size = max_workers * chunk_size

        from multiprocessing import Process
        input_queue = self.manager.Queue(maxsize=queue_max_size)
        output_queue = self.manager.Queue(maxsize=queue_max_size)

        # Initializing processes
        processes = [
            Process(
                target=MultiProcessingQueueExecutor._worker,
                args=(input_queue, output_queue, fn)
            )
            for _ in range(max_workers)
        ]
        for p in processes:
            p.start()

        # Initial filling of the input queue
        all_args = iter(zip(*iterables))
        for _ in range(queue_max_size):
            try:
                input_queue.put(next(all_args))
            except StopIteration:
                input_queue.put(None)
                break

        # Getting and yielding results
        while True:
            result = output_queue.get()
            if result is None:
                max_workers -= 1
                if max_workers == 0:
                    break
            else:
                yield result
                try:
                    input_queue.put(next(all_args))
                except StopIteration:
                    input_queue.put(None)

        # Cleaning
        for p in processes:
            p.join()

    @staticmethod
    def _worker(input_queue: "Queue", output_queue: "Queue", fn) -> None:
        while True:
            args = input_queue.get()
            if args is None:
                break
            result = fn(*args)
            output_queue.put(result)

        input_queue.put(None)
        output_queue.put(None)

    
ParallelExecutorStr = Literal[
    "Sequential",
    "MultiProcessing",
    "FutureProcess",
    "FutureThread",
    "FutureMPI",
    "MultiProcessingQueue",
]

parallel_executor_dict: Final[dict[ParallelExecutorStr, type[ParallelExecutor]]] = {
    "Sequential": SequentialExecutor,
    "MultiProcessing": MultiProcessingPoolExecutor,
    "FutureProcess": FutureProcessExecutor,
    "FutureThread": FutureThreadExecutor,
    "FutureMPI": FutureMPIExecutor,
    "MultiProcessingQueue": MultiProcessingQueueExecutor,
}

class Parallel(AbstractContextManager["Parallel"]):
    """ Unique parallel context
    
    This generates a singleton of a parallel executor
    """
    executor_class: ClassVar[type[ParallelExecutor]] = SequentialExecutor
    max_workers: ClassVar[Optional[int]] = None
    chunk_size: ClassVar[int] = 1
    unordered: ClassVar[bool] = False
    kwargs: ClassVar[dict[str, Any]] = dict()
    __executor: ClassVar[Optional[ParallelExecutor]] = None

    def __init__(self) -> None:
        self.start()

    @staticmethod
    def configure(executor_class: type[ParallelExecutor] | ParallelExecutorStr = "FutureProcess",
                  max_workers: Optional[int] = None,
                  chunk_size: int = 1,
                  unordered: bool = False,
                  **kwargs: Any) -> None:
        """ Configure the parallel executor class and it's init parameters """
        if isinstance(executor_class, str):
            from .utils import to_literal
            executor_class = cast(ParallelExecutorStr, to_literal(ParallelExecutorStr, executor_class))
            executor_class = parallel_executor_dict[executor_class]

        Parallel.executor_class = executor_class
        Parallel.max_workers = max_workers
        Parallel.chunk_size = chunk_size
        Parallel.unordered = unordered
        Parallel.kwargs = kwargs

    @staticmethod
    def start() -> None:
        """ Start a parallel context """
        if Parallel.__executor is None:
            Parallel.__executor = Parallel.executor_class(
                max_workers=Parallel.max_workers,
                chunk_size=Parallel.chunk_size,
                unordered=Parallel.unordered,
                **Parallel.kwargs,
            )

    #def stop(self, *args: Any, **kwargs: Any) -> None:
    @staticmethod
    def shutdown(wait: bool = True) -> None:
        """ Stop a parallel context, eventually wait to finish tasks """
        if Parallel.__executor is not None:
            Parallel.__executor.shutdown(wait)
            Parallel.__executor = None

    @property
    def executor(self) -> ParallelExecutor:
        """ Returns associated parallel executor """
        Parallel.start()
        assert Parallel.__executor is not None
        return Parallel.__executor
    
    def __enter__(self) -> Self:
        """ Entering context """
        self.start()
        return self
    
    def __exit__(self, *args: Any, **kwargs: Any) -> None:
        """ Leaving context """
        self.shutdown()


    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments that configure the parallelization """
        import typing
        from .utils import to_literal
        group = parent_parser.add_argument_group(
            "Parallel computation"
        )
        group.add_argument(
            "--parallel",
            type=lambda s: to_literal(ParallelExecutorStr, s),
            choices=typing.get_args(ParallelExecutorStr),
            nargs='?',
            default="Sequential",
            const="FutureProcess",
            help="Type of parallel executor. FutureProcess is used if no executor is specified."
        )
        group.add_argument(
            "--max_workers",
            type=int,
            default=None,
            help="Number of workers for the parallel computation. Use number of cores if None."
        )
        group.add_argument(
            "--chunk_size",
            type=int,
            default=1,
            help="Number of tasks to pass to each work at once. You should increase this value in order to reduce overhead when there are enough tasks comparing to the number of workers."
        )
        group.add_argument(
            "--unordered",
            action="store_true",
            help="Yield elements when ready, possibly out of order. It may slow down the computation but it is necessary to have a working --lazy in parallel."
        )

    @classmethod
    def from_config(cls: type[Self], config: Namespace, **kwargs: Any) -> "Parallel":
        """ Configure the parallel context from the command-line arguments """
        Parallel.configure(
            executor_class=config.parallel,
            max_workers=config.max_workers,
            chunk_size=config.chunk_size,
            unordered=config.unordered,
        )
        return Parallel()
