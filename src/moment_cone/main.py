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


    parser = argparse.ArgumentParser(
        "Redundant list of inequalities for the moment cone",
        description="""This software computes the moment cone for QMP (Kronecker) and fermion""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )

    Representation.add_arguments(parser)

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
    def compute() -> None:
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

    if config.cprofile is None:
        compute()
    else:
        from .utils import cprofile
        cprofile(compute, file_name=config.cprofile)