__all__ = (

)

from .typing import *

"""
def main(V: Representation,
         filters: Sequence[inequality_filters],
         **kwargs: Any,
        ) -> list[Inequality]:
    return NotImplemented
"""

def main_from_cmd() -> None:
    """ Main entrance from command-line """
    import argparse
    parser = argparse.ArgumentParser(
        "Redundant list of inequalities for the cone",
        description="""
This software compute a redundant list of inequalities for a cone""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )
    parser.add_argument(
        "representation",
        type=str.capitalize,
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

    from .main_steps import ConeStep
    ConeStep.add_arguments(parser)

    # Parsing command-line arguments
    config = parser.parse_args()

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
    inequalities = list(step().all())

    # Printing the inequalities
    print(f"Computed {len(inequalities)} inequalities:")
    for ineq in inequalities:
        print(ineq)

    #Task.print_all()
