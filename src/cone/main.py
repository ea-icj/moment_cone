from .typing import *
from .dimension import Dimension
from .task import Task

from .vector_chooser import *
from .hyperplane_candidates import *
from .list_of_W import *
from .inequality import *
from .tau import *
from .ramification import *
from .Normaliz2python import *
from .dim_stabilizer_K import *

import sys

__all__ = (
    "main",
    "main_from_input",
)

def main(d: Iterable[int] | Dimension,
         tpi_method: str | Method = "probabilistic",
         ram_schub_method: str | Method = "probabilistic",
         ram0_method: str | Method = "probabilistic",
         quiet: bool = False,
         ) -> list[Inequality]:
    """
    Main entrance of the algorithm
    
    Example in 2, 2, 2:
    >>> d = Dimension((2, 2, 2))
    >>> for ineq in main(d, quiet=True):
    ...     print(ineq)
    Inequality(tau  = -2 | 1 0 | 1 0 | 1 0,
               w    =     0 1 | 0 1 | 1 0,
               wtau = -2 | 1 0 | 1 0 | 0 1)

    Example in 3, 3, 3:
    >>> d = Dimension((3, 3, 3))
    >>> for ineq in main(d, quiet=True):
    ...     print(ineq)
    Inequality(tau  = -4 | 2 1 0 | 2 1 0 | 2 1 0,
               w    =     0 2 1 | 1 0 2 | 2 0 1,
               wtau = -4 | 2 0 1 | 1 2 0 | 1 0 2)
    Inequality(tau  = -1 | 0 0 0 | 0 0 0 | 1 1 0,
               w    =     0 1 2 | 0 1 2 | 0 1 2,
               wtau = -1 | 0 0 0 | 0 0 0 | 1 1 0)
    Inequality(tau  = -2 | 1 0 0 | 1 0 0 | 1 1 0,
               w    =     0 1 2 | 1 0 2 | 0 2 1,
               wtau = -2 | 1 0 0 | 0 1 0 | 1 0 1)
    Inequality(tau  = -4 | 2 1 0 | 2 1 0 | 2 1 0,
               w    =     0 1 2 | 1 0 2 | 2 1 0,
               wtau = -4 | 2 1 0 | 1 2 0 | 0 1 2)
    Inequality(tau  = -2 | 1 0 0 | 1 0 0 | 1 1 0,
               w    =     0 1 2 | 2 0 1 | 0 1 2,
               wtau = -2 | 1 0 0 | 0 0 1 | 1 1 0)
    Inequality(tau  = -4 | 2 1 0 | 2 1 0 | 2 1 0,
               w    =     0 1 2 | 0 2 1 | 2 1 0,
               wtau = -4 | 2 1 0 | 2 0 1 | 0 1 2)
    Inequality(tau  = -4 | 2 1 0 | 2 1 0 | 2 1 0,
               w    =     0 2 1 | 0 2 1 | 1 2 0,
               wtau = -4 | 2 0 1 | 2 0 1 | 0 2 1)
    Inequality(tau  = -2 | 1 0 0 | 1 0 0 | 1 1 0,
               w    =     0 1 2 | 0 1 2 | 1 2 0,
               wtau = -2 | 1 0 0 | 1 0 0 | 0 1 1)
    """
    if not isinstance(d, Dimension):
        d = Dimension(d)
    tpi_method = to_method(tpi_method)
    ram_schub_method = to_method(ram_schub_method)
    ram0_method = to_method(ram0_method)

    if not quiet:
        print("d =", d)
        print("tpi_method =", tpi_method)
        print("ram_schub_method =", ram_schub_method)
        print("ram0_method =", ram0_method)
        print()

    # Reset task history
    Task.reset_all()
    Task.quiet = quiet

    ## Checking if the cone has the expected dimension
    with Task("Step 0: checking cone dimension"):
        Ms=Lie_action_as_matrices_V(d)
        MsR=[mat_C_to_R(M) for M in Ms]
        if dim_gen_stab_of_K(MsR)>len(d)-1: # Check that the dim is computed in U_n(C)^s without the isolated S^1
            raise ValueError(
                f'The moment cone has codimension bigger that the length of d. Namely {dim_gen_stab_of_K(MsR)}.'
                'The program does not work in this case.'
            )
        elif not quiet:
            print('The cone has the expected dimension.')
            print()

    ## Generate the list of candidates for tau
    #Candidates_for_tau=find_hyperplanes_mod_sym_dim(d0,d0.dimU) # This is the function for regular ops (todo : include this info in the name) - To be changed.
    with Task('Step 1: looking for a first list of dominant 1-PS whose kernel is supported at hyperplanes of weights'):
        Candidates_for_tau_1ps = find_1PS_mod_sym_dim(d, quiet)

    if not quiet:
        print(len(Candidates_for_tau_1ps), ' dominant 1-PS selected in Step 1')
        for tau in Candidates_for_tau_1ps:
            print(tau)
        print()

    with Task('Step 1.5: unicity of this list'):
        Candidates_for_tau = unique_modulo_symmetry_list_of_tau(Candidates_for_tau_1ps) # todo : inutile car déjà fait dans find_1PS ?

    if not quiet:
        print(len(Candidates_for_tau), ' dominant 1-PS selected in Step 1 After Unicity')
        for tau in Candidates_for_tau:
            print(tau)
        print()
 
    # Filter 1: submodule condition
    with Task('Step 2: checking submodule condition'):
        Candidates_for_tau1=[tau for tau in Candidates_for_tau if tau.is_sub_module]

    if not quiet:
        print(len(Candidates_for_tau1), ' dominant 1-PS selected in Step 2')
        print()

    # Filter 2: stabilizer condition
    with Task('Step 3: stabilizer condition'):
        ### Avec le nouveau dimStab
        Candidates_for_tau2=[]
        for tau in Candidates_for_tau1:
            if  tau.is_dom_reg :
                Candidates_for_tau2.append(tau)
            else: 
                Ms_tau=Lie_action_as_matrices_Vtau(tau,Ms)
                Ms_tauR=[mat_C_to_R(M) for M in Ms_tau]
                if dim_gen_stab_of_K(Ms_tauR)==len(d):
                    Candidates_for_tau2.append(tau)   

    if not quiet:
        print(len(Candidates_for_tau2), ' dominant 1-PS selected in Step 3')
        for tau in Candidates_for_tau2:
            print(tau)
        print()

    ## Generate the list of candidates for the inequalities (pairs tau,w)
    ## Here w has to belong to P^tau and U(w) is tau-isomorphic to V(tau>0)
    with Task('Step 4: computing action of W, in order to get a first list of inequalities containing all the expected ones'):
        Candidates_for_Ineq: list[Inequality] = []
        for tau in Candidates_for_tau2 :
            Lw = ListWs_Mod(tau)
            Candidates_for_Ineq += [Inequality(tau, w) for w in Lw] # Fait-on un dictionnaire tau : Lw ??

    if not quiet:
        print(len(Candidates_for_Ineq), ' inequalities selected in Step 4')
        for ineq in Candidates_for_Ineq :
            print(ineq)
        print()

    # Filter 0: unicity modulo sym(d)
    with Task('Step 5, Reduction modulo symmetries of the dimension vector'):
        Candidates_for_Ineq1=unique_modulo_symmetry_list_of_ineq(Candidates_for_Ineq)

    if not quiet:
        print(len(Candidates_for_Ineq1), ' inequalities selected in Step 5')
        for ineq in Candidates_for_Ineq1 :
            print(ineq.tau)
        print()

    # Filter 1: pi is dominant
    with Task('Step 6: checking dominancy of the map pi'):
        Dominant_Ineq=[ineq for ineq in Candidates_for_Ineq1 if Check_Rank_Tpi(ineq,tpi_method)] 

    if not quiet:
        print(len(Dominant_Ineq), ' inequalities selected in Step 6')
        for ineq in Dominant_Ineq :
            print(ineq.tau)
        print()

    # Filter 2: pi is birational (ramification divisor contracted)
    with Task('Step 7: checking birationality (ramification divisor contracted) of the map pi'):
        Birational_Ineq=[ineq for ineq in Dominant_Ineq if Is_Ram_contracted(ineq,ram_schub_method,ram0_method)]

    if not quiet:
        print(len(Birational_Ineq), ' inequalities selected in Step 7')    
        for ineq in Birational_Ineq :
            print(ineq)
        print()

    #path="/home/bm29130h/Documents/Recherche/Ressources_autres/GDT/Machine Learning/calculs Kron/2 oct/"
    #reference=[Inequality.from_tau(tau) for tau in convert_file_Nout2pyth(path,d0)]
    #dictionary_list_lengths(compare_ineq_candidates_reference_mod_sym_dim(Candidates_for_Ineq1,reference))
    #test fails following fusion 17dec 12h11

    #unique_reference=unique_modulo_symmetry_list_of_ineq(reference)
    #dictionary_list_lengths(compare_tau_candidates_reference_mod_sym_dim(Candidates_for_tau,[ineq.tau for ineq in reference]))

    
    # Disabling this part since the reference files are missing
    # TODO: move to a unittest ?
    if False and d in [Dimension([3,3,3]), Dimension([4,4,4]), Dimension([5,4,4])]:
        print()
        print(80 * "#")
        print("Comparison with values of reference for ", d)
        path="./cone/References/"
        reference=[Inequality.from_tau(tau) for tau in convert_file_Nout2pyth(path,d)]
        diff_with_reference=compare_ineq_candidates_reference_mod_sym_dim( Birational_Ineq , reference)
        print( dictionary_list_lengths(diff_with_reference))
        #FIXME: d[0]-1 has to be adapted when values differ (prod of different d[i]-1)
        print("if perfect: should return 0 for candidates_only and ", d[0]-1, "for reference_only")
    
    if not quiet:
        print()
        print(80 * "#")
        print("Summary of the tasks:")
        Task.print_all(disp_interlude=False)

    return Birational_Ineq

def to_method(method: str) -> Method:
    """ Return full name of a method by looking only to the first letter """
    if method.startswith("p"):
        return "probabilistic"
    elif method.startswith("s"):
        return "symbolic"
    else:
        raise ValueError(f"Unknown method {method}")

def main_from_input() -> list[Inequality]:
    """ Main entrance from user inputs """
    print('This software compute a redundant list of inequalities for the cone Kron(d1, d2, ..., ds) with d1 >= d2 >= d3 ...')

    print()
    d_str = input("Enter the d_i separated by commas: ")
    d = Dimension(map(int, d_str.split(",")))

    print()
    print("For the following steps, you can choose either the probabilistic or symbolic method to be used (please type p or s in each case):")
    tpi_method = to_method(input("\tFor surjectivity of Tpi: "))
    ram_schub_method = to_method(input("\tFor checking if the Bruhat ramification divisors are contracted: "))
    ram0_method = to_method(input("\tFor checking if R0 is contracted: "))

    return main(d,
         tpi_method=tpi_method,
         ram_schub_method=ram_schub_method,
         ram0_method=ram0_method
    )

def main_from_cmd() -> None:
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
    parser.add_argument("--tpi", type=to_method, choices=method_choices, default="p", help="Method for the surjectivity of Tpi (p for probabilistic, s for symbolic)")
    parser.add_argument("--ram_schub", type=to_method, choices=method_choices, default="p", help="Method for checking if the Bruhat ramification divisors are contracted (p for probabilistic, s for symbolic)")
    parser.add_argument("--ram0", type=to_method, choices=method_choices, default="p", help="Method for checking if R0 is contracted (p for probabilistic, s for symbolic)")
    parser.add_argument("--quiet", action="store_true", help="Removing output during computation")

    config = parser.parse_args()

    main(Dimension(config.d),
         tpi_method=config.tpi,
         ram_schub_method=config.ram_schub,
         ram0_method=config.ram0,
         quiet=config.quiet,
    )
