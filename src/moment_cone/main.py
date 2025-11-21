__all__ = (
    'moment_cone',
    'moment_cone_from_cmd',
)

import typing

from .typing import *
from .representation import Representation
from .inequality import Inequality
from .main_steps import Dataset, InequalityFilterStr, default_inequalities_filters
from .main_steps import MomentConeStep

if TYPE_CHECKING:
    from tracemalloc import Snapshot as TraceMallocSnapshot


@overload
def moment_cone(V: Representation,
         filters: Sequence[InequalityFilterStr],
         return_step: Literal[False],
         **options: Any) -> Dataset[Inequality]:
    ...

@overload
def moment_cone(V: Representation,
         filters: Sequence[InequalityFilterStr],
         return_step: Literal[True],
         **options: Any) -> tuple[Dataset[Inequality], MomentConeStep]:
    ...

def moment_cone(V: Representation,
         filters: Sequence[InequalityFilterStr] = default_inequalities_filters,
         return_step: bool = False, # Returns the dataset and the ConeStep instance
         **options: Any
    ) -> Dataset[Inequality] | tuple[Dataset[Inequality], MomentConeStep]:
    """ Main entrance from Python prompt

    For the options, see the description of each step in main_steps.
    """
    cone_step = MomentConeStep(V, filters=filters, **options)
    if return_step:
        return cone_step(), cone_step
    else:
        return cone_step()


def moment_cone_from_cmd() -> None:
    """ Main entrance from command-line """
    import argparse
    from .utils import to_literal
    from .representation import Representation
    from .parallel import Parallel


    parser = argparse.ArgumentParser(
        "Redundant list of inequalities for the moment cone",
        description="""This software computes the moment cone for QMP (Kronecker) and fermion""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )
    parser.add_argument(
        "--check_inequalities",
        action="store_true",
        help="Check computed inequalities with references datas"
    )

    Representation.add_arguments(parser)
    Parallel.add_arguments(parser)

    group = parser.add_argument_group("Development tools")
    group.add_argument(
        "--line_profiler",
        type=str,
        nargs="*",
        default=[],
        help="Profile given function by line",
    )
    group.add_argument(
        "--cprofile",
        type=str,
        default=None,
        help="Profile function calls and output results in given file name (pstats and kcachegrind format)",
    )
    group.add_argument(
        "--tracemalloc",
        action="store_true",
        help="Use tracemalloc library to trace memory allocation during the computation"
    )
    group.add_argument(
        "--tm_frame",
        type=int,
        default=1,
        help="Number of frames to consider when using tracemalloc"
    )
    group.add_argument(
        "--tm_top",
        type=int,
        default=10,
        help="Number of lines to display when using tracemalloc",
    )
    group.add_argument(
        "--logging_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        type=str.upper,
        default="INFO",
        help="Level of the displayed logging messages",
    )

    from .main_steps import MomentConeStep
    MomentConeStep.add_arguments(parser)

    # Parsing command-line arguments
    config = parser.parse_args()

    # Seed
    from .utils import manual_seed
    config.seed = manual_seed(config.seed)

    # Parallel context
    Parallel.from_config(config)

    # Configuring the logging level
    from .utils import getLogger
    import logging
    getLogger().setLevel(getattr(logging, config.logging_level))

    # Displaying configuration
    if not config.quiet:
        print("Configuration:")
        for k, v in vars(config).items():
            print(f"\t{k}: {v}")
        print()

    # Creating the representation
    V = Representation.from_config(config)
    
    # Creating the overall cone computational step
    step = MomentConeStep.from_config(V, config)

    # Reset task history
    from .task import Task
    Task.reset_all()

    # Computing the cone
    def compute() -> list[Inequality]:
        if len(config.line_profiler) > 0:
            from .utils import line_profiler
            inequalities, lp = line_profiler(
                config.line_profiler,
                lambda: list(step())
            )

            print("\nLine profiling results:")
            lp.print_stats()
        else:
            inequalities = list(step())
        return inequalities

    if config.tracemalloc:
        import tracemalloc
        tracemalloc.start(config.tm_frame)

    if config.cprofile is None:
        inequalities = compute()
    else:
        from .utils import cprofile
        inequalities, stats = cprofile(compute, file_name=config.cprofile)

    if config.tracemalloc:
        import tracemalloc
        snapshot = tracemalloc.take_snapshot()
        tracemalloc.stop()
        display_top(snapshot, limit=config.tm_top)



    # Checking inequalities
    if config.check_inequalities:
        from .reference_datas import compare_to_reference
        print()
        try:
            cmp = compare_to_reference(inequalities, V, source="results")
        except KeyError:
            print(f"There are no reference inequalities for {V}")
        else:
            print("Comparison with the reference inequalities:")
            print(cmp)
            if not cmp:
                print("Only in results:")
                for ineq in cmp.only1:
                    print(f"\t{ineq}")
                print("Only in reference:")
                for ineq in cmp.only2:
                    print(f"\t{ineq}")
                

# From the official documentation: https://docs.python.org/3/library/tracemalloc.html#pretty-top
def display_top(snapshot: "TraceMallocSnapshot", key_type: str = 'lineno', limit: int = 10) -> None:
    import tracemalloc
    import linecache

    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB"
              % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))