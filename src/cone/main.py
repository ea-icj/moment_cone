from cone.typing import *
from cone.dimension import Dimension

def to_method(method: str) -> Method:
    if method.startswith("p"):
        return "probabilistic"
    elif method.startswith("s"):
        return "symbolic"
    else:
        raise ValueError(f"Unknown method {method}")

def main(d: Iterable[int] | Dimension,
         stabilizer_method: str | Method,
         tpi_method: str | Method,
         ram_schub_method: str | Method,
         ram0_method: str | Method,
         ):
    """ Main entrance of the algorithm """
    if not isinstance(d, Dimension):
        d = Dimension(d)
    stabilizer_method = to_method(stabilizer_method)
    tpi_method = to_method(tpi_method)
    ram_schub_method = to_method(ram_schub_method)
    ram0_method = to_method(ram0_method)

    print("d =", d)
    print("stabilizer_method =", stabilizer_method)
    print("tpi_method =", tpi_method)
    print("ram_schub_method =", ram_schub_method)
    print("ram0_method =", ram0_method)

    pass # TODO

def main_from_input():
    """ Main entrance from user inputs """
    print('This software compute a redundant list of inequalities for the cone Kron(d1, d2, ..., ds) with d1 >= d2 >= d3 ...')

    print()
    d_str = input("Enter the d_i separated by commas: ")
    d = Dimension(map(int, d_str.split(",")))

    print()
    print("For the following steps, you can choose either the probabilistic or symbolic method to be used (please type p or s in each case):")
    stabilizer_method = to_method(input("\tFor stabilizer in K: "))
    tpi_method = to_method(input("\tFor surjectivity of Tpi: "))
    ram_schub_method = to_method(input("\tFor checking if the Bruhat ramification divisors are contracted: "))
    ram0_method = to_method(input("\tFor checking if R0 is contracted: "))

    main(d,
         stabilizer_method=stabilizer_method,
         tpi_method=tpi_method,
         ram_schub_method=ram_schub_method,
         ram0_method=ram0_method
    )

def main_from_cmd():
    """ Main entrance from command-line """
    import argparse
    parser = argparse.ArgumentParser(
        "Redundant list of inequalities for the cone Kron",
        description="""
This software compute a redundant list of inequalities for the cone Kron(d1, d2, ..., ds) with d1 >= d2 >= d3 ...
""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,    
    )

    parser.add_argument("d", type=int, nargs="+", help="The dimensions d_i")

    method_choices = ["p", "s", "probabilistic", "symbolic"]
    parser.add_argument("--stabilizer", type=to_method, choices=method_choices, default="p", help="Method for the stabilizer in K (p for probabilistic, s for symbolic)")
    parser.add_argument("--tpi", type=to_method, choices=method_choices, default="p", help="Method for the surjectivity of Tpi (p for probabilistic, s for symbolic)")
    parser.add_argument("--ram_schub", type=to_method, choices=method_choices, default="p", help="Method for checking if the Bruhat ramification divisors are contracted (p for probabilistic, s for symbolic)")
    parser.add_argument("--ram0", type=to_method, choices=method_choices, default="p", help="Method for checking if R0 is contracted (p for probabilistic, s for symbolic)")

    config = parser.parse_args()

    main(Dimension(config.d),
         stabilizer_method=config.stabilizer,
         tpi_method=config.tpi,
         ram_schub_method=config.ram_schub,
         ram0_method=config.ram0,
    )

if __name__ == "__main__":
    main_from_cmd()
