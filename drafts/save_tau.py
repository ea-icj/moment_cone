#!/usr/bin/env python3
"""
Backup the tau generated after StabilizerCondition step
"""
def main_from_cmd() -> None:
    import argparse
    from moment_cone import Representation
    from moment_cone import Step, TauCandidatesStep, SubModuleConditionStep, StabilizerConditionStep
    from moment_cone.task import Task
    from moment_cone.parallel import Parallel

    parser = argparse.ArgumentParser(
        "Backup the tau generated after StabilizerCondition step",
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
        help="File name where to store the results",
    )

    Representation.add_arguments(parser)
    Parallel.add_arguments(parser)
    TauCandidatesStep.add_arguments(parser)
    SubModuleConditionStep.add_arguments(parser)
    StabilizerConditionStep.add_arguments(parser)

    # Parsing command-line arguments
    config = parser.parse_args()

    # Seed
    from moment_cone.utils import manual_seed
    config.seed = manual_seed(config.seed)

    # Parallel context
    Parallel.from_config(config)

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
    def create_step(step_type: type[TStep]) -> TStep:
        """ Create and configure a new step """
        dataset_type = LazyDataset if config.lazy else ListDataset
        return step_type.from_config(V, config, dataset_type=dataset_type)
    
    # Tau candidates
    tau_candidates_step = create_step(TauCandidatesStep)
    with Task(tau_candidates_step.name) as task:
        tau_candidates = tau_candidates_step()
        task.log(f"tau_candidates: {tau_candidates}", indent=1)


    # Filters candidate tau
    for tau_filter_type in SubModuleConditionStep, StabilizerConditionStep:
        tau_filter_step = create_step(tau_filter_type)
        with Task(tau_filter_step.name) as task:
            tau_candidates = tau_filter_step(tau_candidates)
            task.log(f"tau_candidates: {tau_candidates}", indent=1)
    
    # CSV export
    if config.file_name is None:
        config.file_name = f"tau_{V}"
    
    import csv
    from tqdm.auto import tqdm
    with Task("Export") as task:
        fh = open(f"{config.file_name}.csv", "w")
        csvw = csv.writer(fh, quoting=csv.QUOTE_NONNUMERIC)
        for tau in tqdm(tau_candidates.pending(), desc="Export", unit="tau"):
            csvw.writerow(tau.flattened)
            fh.flush()
        fh.close()

if __name__ == "__main__":
    main_from_cmd()
