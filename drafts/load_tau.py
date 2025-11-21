#!/usr/bin/env python3
"""
Reload saved tau and resume computation after StabilizerCondition step
"""

def main_from_cmd() -> None:
    import argparse
    import typing
    from typing import Iterable, Any

    from moment_cone import Representation
    from moment_cone import Step, InequalityCandidatesStep, ExportStep
    from moment_cone.main_steps import InequalityFilterStr, inequalities_filter_dict, default_inequalities_filters
    from moment_cone.utils import to_literal
    from moment_cone.parallel import Parallel

    parser = argparse.ArgumentParser(
        "Reload saved tau and resume computation after StabilizerCondition step",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable informations during computation",
    )
    parser.add_argument(
        "--lazy",
        action="store_true",
        help="Compute lazilly the inequalities (without storing intermediate results)",
    )
    parser.add_argument(
        "--file_name",
        type=str,
        default=None,
        help="File name from where to load the tau",
    )
    parser.add_argument( # FIXME: Move in step
        "--filters",
        type=lambda s: to_literal(InequalityFilterStr, s),
        nargs='*',
        choices=typing.get_args(InequalityFilterStr),
        default=default_inequalities_filters,
        help="Sequence of filters applied to the inequalities",
    )
    parser.add_argument(
        "--check_inequalities",
        action="store_true",
        help="Check computed inequalities with references datas"
    )
    parser.add_argument(
        "--logging_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        type=str.upper,
        default="INFO",
        help="Level of the displayed logging messages",
    )


    Representation.add_arguments(parser)
    Parallel.add_arguments(parser)
    InequalityCandidatesStep.add_arguments(parser)

    # Adding command-line options from inequality filtering steps
    for name, cls in inequalities_filter_dict.items():
        cls.add_arguments(parser)
    ExportStep.add_arguments(parser)

    # Parsing command-line arguments
    config = parser.parse_args()

    # Seed
    from moment_cone.utils import manual_seed
    config.seed = manual_seed(config.seed)

    # Parallel context
    Parallel.from_config(config)
    
    # Configuring the logging level
    from moment_cone.utils import getLogger
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

    # Helper to create a step
    from typing import TypeVar
    from moment_cone.main_steps import ListDataset, LazyDataset
    TStep = TypeVar("TStep", bound=Step)
    def create_step(step_type: type[TStep], **kwargs: Any) -> TStep:
        """ Create and configure a new step """
        dataset_type = LazyDataset if config.lazy else ListDataset
        return step_type.from_config(V, config, dataset_type=dataset_type, **kwargs)

    # Initial dataset with tau from loaded file
    from moment_cone import Tau
    if config.file_name is None:
        config.file_name = f"tau_{V}.csv"

    # Loading all tau
    import csv
    fh = open(config.file_name, "r")
    csvr = csv.reader(fh, quoting=csv.QUOTE_NONNUMERIC)
    all_tau = [Tau.from_flatten(map(int, row), V.G) for row in csvr]
    fh.close()

    from moment_cone.task import Task
    from moment_cone import Inequality, Dataset
    from tqdm.auto import tqdm

    with Task("Import") as task:
        dataset_type = LazyDataset[Tau] if config.lazy else ListDataset[Tau]
        tau_candidates = dataset_type.from_separate(
            tqdm(all_tau,  desc="Import", unit="tau")
        )

    # Transform tau to inequality
    ineq_candidates: Dataset[Inequality]
    ineq_candidates_step = create_step(InequalityCandidatesStep)
    with Task(ineq_candidates_step.name) as task:
        ineq_candidates = ineq_candidates_step(tau_candidates)
        task.log(f"ineq_candidates: {ineq_candidates}", indent=1)
    
    # Pre-computation of Representation.TPi 3D matrix if necessary
    if 'PiDominancy' in config.filters or 'Birationality' in config.filters:
        from moment_cone.main_steps import TPiPreComputationStep
        TPi_step = create_step(TPiPreComputationStep)
        with Task(TPi_step.name):
            ineq_candidates = TPi_step(ineq_candidates)

    # Filters candidate inequalities
    from pathlib import Path
    import csv
    fh_mult = open(str(Path(config.file_name).stem) + "_mult.csv", "w")
    csvw = csv.writer(fh_mult, quoting=csv.QUOTE_NONNUMERIC)
    
    for name in config.filters:
        ineq_filter_type = inequalities_filter_dict[name]
        if name == "Birationality":
            def mult_callback(ineq: Inequality, mult: int):
                csvw.writerow((mult,) + tuple(ineq.wtau.flattened))
                fh_mult.flush()
            kwargs = dict(mult_callback=mult_callback)
        else:
            kwargs = dict()
        ineq_filter_step = create_step(ineq_filter_type, **kwargs)
        with Task(ineq_filter_step.name) as task:
            ineq_candidates = ineq_filter_step(ineq_candidates)
            task.log(f"ineq_candidates: {ineq_candidates}", indent=1)
    
    # Exporting inequalities
    export_step = create_step(ExportStep)
    with Task(export_step.name):
        ineq_candidates = export_step(ineq_candidates)    

    fh_mult.close()

if __name__ == "__main__":
    main_from_cmd()
