import contextlib
import itertools
import functools
import time
import signal
from psutil import Process

from .typing import *

__all__ = (
    "Task",
    "TimeOutException",
    "timeout",
)

class Task(contextlib.AbstractContextManager):
    """
    Context manager to measure and log task durations

    Remark that cpu time measurement may be inconsistent if processes
    are created or destructed during the usage of this class.

    Example:
    with Task("Computing stuff"):
        # do things
        ...
        ...

    Task.print_all()
    """
    name: str
    perf_counter: tuple[Optional[int], Optional[int]]
    process_time: tuple[Optional[int], Optional[int]]

    process: Process = Process()
    all_tasks: list["Task"] = [] # All created tasks (static)
    all_start: tuple[int, int] = (0, 0)
    quiet: bool = False


    @classmethod
    def current_wall_time(cls) -> int:
        """ Wall time in ns """
        return time.perf_counter_ns()
    
    @classmethod
    def current_process_time(cls) -> int:
        """ Process (including childrens) CPU time (user + system) """
        cpu_times = sum(cpu_times.user + cpu_times.system for p in cls.process.children(recursive=True) for cpu_times in (p.cpu_times(),))
        return round(1e9 * cpu_times) + time.process_time_ns()
    
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
            stop = (cls.current_wall_time(), cls.current_process_time())
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
        cls.all_start = (cls.current_wall_time(), cls.current_process_time())

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

    def __init__(self, name: str, auto_start: bool = False):
        """ Manual construction """
        self.name = name
        self.perf_counter = None, None
        self.process_time = None, None
        self.all_tasks.append(self)
        if auto_start:
            self.start()

    def __enter__(self):
        """ Entering context """
        self.start()
        if not self.quiet:
            print(self, end='\r')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """ Leaving context """
        self.stop()
        if not self.quiet:
            print(self)

    def start(self) -> None:
        assert Task.is_clear(self)
        self.perf_counter = (self.current_wall_time(), None)
        self.process_time = (self.current_process_time(), None)

    def stop(self) -> None:
        assert Task.is_running(self)
        self.perf_counter = (self.perf_counter[0], self.current_wall_time())
        self.process_time = (self.process_time[0], self.current_process_time())

    @property
    def duration(self) -> tuple[int, int]:
        if Task.is_finished(self):
            return (
                self.perf_counter[1] - self.perf_counter[0],
                self.process_time[1] - self.process_time[0]
            )
        elif Task.is_running(self):
            return (
                self.current_wall_time() - self.perf_counter[0],
                self.current_process_time() - self.process_time[0]
            )
        else:
            return (0, 0)

    def __repr__(self) -> str:
        if Task.is_clear(self):
            status = "Waiting"
        elif Task.is_running(self):
            status = "..."
        else:
            status = "Done"

        pc, pt = self.duration
        return f"{self.name}: {status} ({Task.format_wall_cpu(self.duration)})"
    

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
class timeout:
    """
    Decorator and context manager to limit wall execution time of a code
    
    Example of usage as a decorator:

    ```
    @timeout(10)
    def compute(a, b, c):
        print(a, b, c)
        ... # do some stuff
        return a * b + c
    result = compute(1, 2, 3) # result is TimeOutException if task didn't finished
    ```

    Example of usage as a context manager:

    ```
    a, b, c = 1, 2, 3
    with timeout(10):
        print(a, b, c)
        ... # do some stuff
        result = a * b + c
    # result isn't defined at all if time is out
    ```

    Example of usage as a context manager with raised exception:

    ```
    a, b, c = 1, 2, 3
    try:
        with timeout(10)
            print(a, b, c)
            ... # do some stuff
            result = a * b + c
    except TimeOutException:
        # Some something when task didn't finished
    ```
    """
    t: int
    no_raise: bool

    def __init__(self, t: int, no_raise: bool = True):
        """ Initialization of the timeout.
        
        Computational time in seconds.
        no_raise to False if you want to capture the exception (usefull in a context)
        """
        self.t = t
        self.no_raise = no_raise

    def handler(self, signum, frame) -> None:
        """ Signal handler called when time is out """
        raise TimeOutException("Time is out!")
    
    def __call__(self, fn: Callable) -> Callable:
        """
        Decorator interface
        
        Return value is a TimeOutException if execution didn't finished
        before the time is out.
        """
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            result = TimeOutException("Time is out!")
            with self:
                result = fn(*args, **kwargs)
            return result
        
        return wrapper
    
    def __enter__(self):
        """ Beginning of the context """
        signal.signal(signal.SIGALRM, self.handler)
        signal.alarm(self.t)

    def __exit__(self, exc_type, exc_value, traceback):
        """ End of the context """
        signal.alarm(0)
        return self.no_raise and exc_type == TimeOutException
