__all__ = (
    'cone',
    'cone_from_cmd',
)

import typing

from .typing import *
from .representation import Representation
from .inequality import Inequality
from .main_steps import Dataset, InequalityFilterStr, default_inequalities_filters
from .main_steps import ConeStep

@overload
def cone(V: Representation,
         filters: Sequence[InequalityFilterStr],
         return_step: Literal[False],
         **options: Any) -> Dataset[Inequality]:
    ...

@overload
def cone(V: Representation,
         filters: Sequence[InequalityFilterStr],
         return_step: Literal[True],
         **options: Any) -> tuple[Dataset[Inequality], ConeStep]:
    ...

def cone(V: Representation,
         filters: Sequence[InequalityFilterStr] = default_inequalities_filters,
         return_step: bool = False, # Returns the dataset and the ConeStep instance
         **options: Any
    ) -> Dataset[Inequality] | tuple[Dataset[Inequality], ConeStep]:
    """ Main entrance from Python prompt

    For the options, see the description of each step in main_steps.
    """
    cone_step = ConeStep(V, filters=filters, **options)
    if return_step:
        return cone_step(), cone_step
    else:
        return cone_step()


def cone_from_cmd() -> None:
    """ Main entrance from command-line """
    import argparse
    from .utils import to_literal

    parser = argparse.ArgumentParser(
        "Redundant list of inequalities for the cone",
        description="""This software compute a irredundant list of inequalities for a cone""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )
    parser.add_argument(
        "representation",
        type=lambda s: to_literal(Literal["Kronecker", "Fermion", "Boson"], s),
        choices=("Kronecker", "Fermion", "Boson"),
        help="Representation type",
    )
    parser.add_argument(
        "N",
        type=int,
        nargs='+',
        help="Dimensions of the linear groups",
    )
    parser.add_argument(
        "--particle_cnt",
        type=int,
        default=None,
        help="Number of particles in Fermion and Boson representation",
    )
    parser.add_argument(
        "--line_profiler",
        type=str,
        nargs="*",
        default=[],
        help="Profile given function by line",
    )
    parser.add_argument(
        "--cprofile",
        type=str,
        default=None,
        help="Profile function calls and output results in given file name (pstats and kcachegrind format)",
    )
    parser.add_argument(
        "--logging_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        type=str.upper,
        default="INFO",
        help="Level of the displayed logging messages",
    )

    from .main_steps import ConeStep
    ConeStep.add_arguments(parser)

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
    from .linear_group import LinearGroup
    from .representation import Representation
    G = LinearGroup(config.N)
    V: Representation
    match config.representation.lower():
        case "kronecker":
            from .representation import KroneckerRepresentation
            V = KroneckerRepresentation(G)
        case "fermion":
            assert config.particle_cnt is not None, "particle_cnt is mandatory for Fermion representation"
            from .representation import FermionRepresentation
            V = FermionRepresentation(G, particle_cnt=config.particle_cnt)
        case "boson":
            assert config.particle_cnt is not None, "particle_cnt is mandatory for Boson representation"
            from .representation import BosonRepresentation
            V = BosonRepresentation(G, particle_cnt=config.particle_cnt)
        case _:
            raise ValueError(f"Invalid representation name {config.representation}")
    
    # Creating the overall cone computational step
    step = ConeStep.from_config(V, config)

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