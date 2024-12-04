"""
Basic tools to manipulate polynomial rings with variables
"""
from sage.all import Ring # type: ignore
from .typing import *
from .weight import Weight

__all__ = (
    'VariableName',
    'PolynomialRingForWeights',
    'Ring',
)

VariableName = str | Weight
Variable = Any
RingGens = dict[str, Variable]

def variable_name(name_or_weight: VariableName, seed: str = "v") -> str:
    """ Accepts a variable name or a sequence of integer (typically a weight) and returns the corresponding variable name """
    if isinstance(name_or_weight, str):
        return name_or_weight
    else:
        return seed + "_" + "_".join(map(str, name_or_weight))

def variable(ring_or_gens: Ring | RingGens, name_or_weight: VariableName, seed: str = "v") -> Variable:
    """ Get variable of a ring from it's name or sequence of coefficients (typically a weight) """
    return variables(ring_or_gens, (name_or_weight,), seed)[0]
    
def variables(ring_or_gens: Ring | RingGens, name_or_weight: Iterable[VariableName], seed: str = "v") -> tuple[Variable, ...]:
    """ Get multiple variables of a ring from it's name or sequence of coefficients (typically a weight) """
    if not isinstance(ring_or_gens, dict):
        ring_or_gens = ring_or_gens.gens_dict()
    return tuple(ring_or_gens[variable_name(nc, seed)] for nc in name_or_weight)


class PolynomialRingForWeights:
    """
    Aggregate a Sage PolynomialRing with dedicated function to manage variables associated to weights.

    Examples:
    >>> from sage.all import QQ, I
    >>> QZ = PolynomialRingForWeights(QQ, "z")
    >>> QZ
    Univariate Polynomial Ring in z over Rational Field
    >>> QZ.variable('z')
    z

    >>> from cone import Weight, Dimension
    >>> d = Dimension((4, 3, 2))
    >>> weights = list(Weight.all(d))[5:12]
    >>> QV = PolynomialRingForWeights(QQ, weights=weights)
    >>> QV
    Multivariate Polynomial Ring in v_0_2_1, v_1_0_0, v_1_0_1, v_1_1_0, v_1_1_1, v_1_2_0, v_1_2_1 over Rational Field
    >>> QV.variable(weights[2])
    v_1_0_1

    >>> QIV = PolynomialRingForWeights(QQ[I], weights=weights, seed=('vr', 'vi'))
    >>> QIV
    Multivariate Polynomial Ring in vr_0_2_1, vi_0_2_1, vr_1_0_0, vi_1_0_0, vr_1_0_1, vi_1_0_1, vr_1_1_0, vi_1_1_0, vr_1_1_1, vi_1_1_1, vr_1_2_0, vi_1_2_0, vr_1_2_1, vi_1_2_1 over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
    >>> QIV.variable(weights[2])
    (vr_1_0_1, vi_1_0_1)
    """
    sage_ring: Ring
    ring_gens: RingGens
    seed: tuple[str, ...]

    def __init__(self, 
                 base_ring: Ring,
                 names: str | Iterable[str] = (), # Fixed variable names
                 weights: Iterable[Weight] = (), # Generates variables for each weight
                 seed: str | Iterable[str] = "v", # Seed(s) of the variable associated to each weight

                 ):
        if isinstance(names, str):
            names = (names,)
        
        if isinstance(seed, str):
            seed = (seed,)
        else:
            seed = tuple(seed)
        assert len(seed) > 0

        from itertools import chain
        variables_names = list(names) + [
            variable_name(chi, seed=s)
            for chi in weights
            for s in seed
        ]

        from sage.all import PolynomialRing
        self.sage_ring = PolynomialRing(base_ring, variables_names)
        self.ring_gens = self.sage_ring.gens_dict()
        self.seed = seed

    def variable(self, name_or_weight: VariableName) -> Variable:
        if isinstance(name_or_weight, str):
            return variable(self.ring_gens, name_or_weight)
        else:
            variables = tuple(
                variable(self.ring_gens, name_or_weight, seed=s)
                for s in self.seed
            )
            if len(variables) == 1:
                return variables[0]
            else:
                return variables

    def __repr__(self) -> str:
        return repr(self.sage_ring)
    
    def __call__(self, *args, **kwargs) -> Any:
        return self.sage_ring(*args, **kwargs)
    
    def base_ring(self) -> Ring:
        return self.sage_ring.base_ring()
    
    def fraction_field(self) -> Ring:
        return self.sage_ring.fraction_field()

