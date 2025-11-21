#!/usr/bin/env python3
from moment_cone import TauCandidatesStep

def fancy_mem(memory: float) -> str:
    for unit in ('B', 'kiB', 'MiB', 'GiB'):
        if memory <= 1024: break
        memory /= 1024
    return f"{memory:.4g}{unit}"
    

def bench(step: TauCandidatesStep, display_step: int = 100) -> tuple[float, float]:
    import psutil
    import os
    import time
    from moment_cone.tau import find_1PS

    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info()[0]
    initial_time = time.perf_counter()
    max_memory: int = 0

    def get_mem() -> int:
        return process.memory_info()[0] - initial_memory
    
    def get_time() -> float:
        return time.perf_counter() - initial_time
    
    with step._tqdm(step.apply().all(), unit="tau", leave=True) as pb:
        for i, _ in enumerate(pb):
            if i % display_step == 0:
                mem = get_mem()
                t = get_time()
                max_memory = max(max_memory, mem)
                pb.set_postfix(
                    memory=fancy_mem(mem),
                    max_memory=fancy_mem(max_memory),
                )

    return get_time(), max_memory


def main_from_cmd() -> None:
    import argparse
    from moment_cone import Representation
    from moment_cone.main_steps import TauCandidatesStep, LazyDataset

    parser = argparse.ArgumentParser(
        "Memory benchmark of Tau candidates generator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )

    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable informations during computation",
    )

    parser.add_argument(
        "--display_step",
        type=int,
        default=100,
        help="Measure and display metrics each given number of generated tau"
    )
    
    Representation.add_arguments(parser)
    TauCandidatesStep.add_arguments(parser)

    # Parsing command-line arguments
    config = parser.parse_args()

    # Displaying configuration
    if not config.quiet:
        print("Configuration:")
        for k, v in vars(config).items():
            print(f"\t{k}: {v}")
        print()

    # Creating the representation
    V = Representation.from_config(config)
    
    # Creating the overall cone computational step
    step = TauCandidatesStep.from_config(V, config, dataset_type=LazyDataset)

    wall_time, memory = bench(step, display_step=config.display_step)
    print()
    print(f"Wall time: {wall_time:.4g}s")
    print("Memory usage:", fancy_mem(memory))


if __name__ == "__main__":
    main_from_cmd()


