import contextlib
import itertools
import time
from contextlib import contextmanager
import logging

from .typing import *
from .utils import getLogger

__all__ = (
    "Task",
    "TimeOutException",
    "timeout",
)

class Task(contextlib.AbstractContextManager["Task"]):
    """
    Context manager to measure and log task durations

    Example:

    >>> with Task("Computing stuff"):
    ...    # do things
    ...    ...
    ...    ...

    >>> Task.print_all()
    """
    name: str
    perf_counter: tuple[Optional[int], Optional[int]]
    process_time: tuple[Optional[int], Optional[int]]
    level: int

    all_tasks: ClassVar[list["Task"]] = [] # All created tasks (static)
    all_start: ClassVar[tuple[int, int]] = (time.perf_counter_ns(), time.process_time_ns())
    quiet: ClassVar[bool] = False

    @staticmethod
    def is_clear(task: "Task") -> TypeGuard["ClearTask"]:
        """ Is a task clear """
        pc1, pc2 = task.perf_counter
        return pc1 is None and pc2 is None
    
    @staticmethod
    def is_running(task: "Task") -> TypeGuard["RunningTask"]:
        """ Is a task running """
        pc1, pc2 = task.perf_counter
        return pc1 is not None and pc2 is None
        
    @staticmethod
    def is_finished(task: "Task") -> TypeGuard["FinishedTask"]:
        """ Is a task finished """
        pc1, pc2 = task.perf_counter
        return pc1 is not None and pc2 is not None

    @classmethod
    def interlude(cls, t1: Optional["Task"], t2: Optional["Task"]) -> Optional[tuple[int, int]]:
        """ Compute interlude duration between two tasks """
        if t1 is not None and not Task.is_finished(t1):
            return None

        start: tuple[int, int]
        if t1 is None:
            start = cls.all_start
        else:
            start = (t1.perf_counter[1], t1.process_time[1])

        stop: tuple[int, int]
        if t2 is None or not (Task.is_running(t2) or Task.is_finished(t2)):
            stop = (time.perf_counter_ns(), time.process_time_ns())
        else:
            stop = (t2.perf_counter[0], t2.process_time[0])

        return (stop[0] - start[0], stop[1] - start[1])

    @staticmethod
    def format_duration(dt: int) -> str:
        """ Format a duration expressed in ns """
        return f"{dt * 1e-6:.3f}ms"
    
    @staticmethod
    def format_wall_cpu(duration: tuple[int, int]) -> str:
        """ Format a pair of duration for Wall time and CPU time """
        wall, cpu = duration
        format = f"Wall: {Task.format_duration(wall)}, CPU: {Task.format_duration(cpu)}"
        if wall > 0:
            percent = round(100 * cpu / wall)
            format += f" ({percent}%)"
        return format

    @classmethod
    def reset_all(cls) -> None:
        """ Clear the task list """
        cls.all_tasks.clear()
        cls.all_start = (time.perf_counter_ns(), time.process_time_ns())

    @classmethod
    def print_all(cls, disp_interlude: bool = True) -> None:
        """ Display all tasks """
        total_tasks = 0, 0
        total_interludes = 0, 0

        for t1, t2 in itertools.pairwise([None,] + cls.all_tasks + [None,]):
            if t1 is not None:
                print(t1)
                duration = t1.duration
                total_tasks = total_tasks[0] + duration[0], total_tasks[1] + duration[1]

            interlude = Task.interlude(t1, t2)
            if interlude is not None:
                if disp_interlude:
                    print(f"\tinterlude ({Task.format_wall_cpu(interlude)})")
                total_interludes = total_interludes[0] + interlude[0], total_interludes[1] + interlude[1]

        print()
        print(f"Total of {len(cls.all_tasks)} tasks: {Task.format_wall_cpu(total_tasks)}")
        print(f"Total of interludes: {Task.format_wall_cpu(total_interludes)}")

    def __init__(self, name: str, auto_start: bool = False, level: Optional[int] = None):
        """ Manual construction """
        self.name = name
        self.perf_counter = None, None
        self.process_time = None, None
        if level is None:
            level = sum(1 for t in self.all_tasks if Task.is_running(t))
        self.level = level

        self.all_tasks.append(self)
        if auto_start:
            self.start()

    def __enter__(self) -> Self:
        """ Entering context """
        self.start()
        if not self.quiet:
            self.self_log(format="{name}: {status}")
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        """ Leaving context """
        self.stop()
        self.self_log(format="{status} ({duration})", indent=1)

    def start(self) -> None:
        assert Task.is_clear(self)
        self.perf_counter = (time.perf_counter_ns(), None)
        self.process_time = (time.process_time_ns(), None)

    def stop(self) -> None:
        assert Task.is_running(self)
        self.perf_counter = (self.perf_counter[0], time.perf_counter_ns())
        self.process_time = (self.process_time[0], time.process_time_ns())

    @property
    def duration(self) -> tuple[int, int]:
        if Task.is_finished(self):
            return (
                self.perf_counter[1] - self.perf_counter[0],
                self.process_time[1] - self.process_time[0]
            )
        elif Task.is_running(self):
            return (
                time.perf_counter_ns() - self.perf_counter[0],
                time.process_time_ns() - self.process_time[0]
            )
        else:
            return (0, 0)

    def status_str(self) -> str:
        if Task.is_clear(self):
            return "Waiting"
        elif Task.is_running(self):
            return "..."
        else:
            return "Done"        
    
    def __repr__(self) -> str:
        status = self.status_str()
        return f"{self.name}: {status} ({Task.format_wall_cpu(self.duration)})"
    
    def self_log(self, level: int = logging.INFO, format: str = "{name}: {status} ({duration})", indent: int = 0) -> None:
        name = self.name
        status = self.status_str()
        duration = Task.format_wall_cpu(self.duration)
        msg = format.format(name=name, status=status, duration=duration)
        self.log(msg, level, indent)

    def log(self, msg: str, level: int = logging.INFO, indent: int = 0) -> None:
        if not self.quiet:
            logger = getLogger(self.name, indentation_level=self.level + indent)
            logger.log(level, msg)



# Sub-classes to make type check working
class ClearTask(Task):
    perf_counter: tuple[Optional[int], Optional[int]]
    process_time: tuple[Optional[int], Optional[int]]

class RunningTask(Task):
    perf_counter: tuple[int, Optional[int]]
    process_time: tuple[int, Optional[int]]

class FinishedTask(Task):
    perf_counter: tuple[int, int]
    process_time: tuple[int, int]


class TimeOutException(Exception):
    """ Exception raised when a task reach the assigned wall time """
    pass


@contextmanager
def timeout(t: float, no_raise: bool = True) -> Generator[None]:
    """
    Decorator and context manager to limit wall execution time of a code
    
    Negative or zero timeout disable the execution time.

    Example of usage as a decorator:

    >>> @timeout(10)
    ... def compute(a, b, c):
    ...     print(a, b, c)
    ...     pass # Do some computations
    ...     return a * b + c
    ...
    >>> result = compute(1, 2, 3) # result is TimeOutException if task didn't finished
    1 2 3

    Example of usage as a context manager:

    >>> a, b, c = 1, 2, 3
    >>> with timeout(10):
    ...    print(a, b, c)
    ...    pass # do some stuff
    ...    result = a * b + c # result isn't defined at all if time is out
    1 2 3

    Example of usage as a context manager with raised exception:
    
    >>> a, b, c = 1, 2, 3
    >>> try:
    ...    with timeout(10, no_raise=False):
    ...        print(a, b, c)
    ...        pass # do some stuff
    ...        result = a * b + c
    ... except TimeOutException:
    ...     pass # Some something when task didn't finished
    1 2 3
    """
    if t <= 0:
        yield
    else:
        from cysignals.alarm import alarm, AlarmInterrupt, cancel_alarm # type: ignore
        try:
            alarm(t)
            yield
        except AlarmInterrupt:
            if not no_raise:
                raise TimeOutException("Time is out!")
        finally:
            cancel_alarm()
