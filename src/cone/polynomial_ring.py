"""
Basic tools to manipulate polynomial rings with variables
"""
from sage.all import Ring # type: ignore
from .typing import *

__all__ = (
    'variable_name',
    'variable',
    'VariableName',
)

VariableName = str | Iterable[int]
Variable = Any

def variable_name(name_or_coeffs: VariableName, seed: str = "v") -> str:
    """ Accepts a variable name or a sequence of integer (typically a weight) and returns the corresponding variable name """
    if isinstance(name_or_coeffs, str):
        return name_or_coeffs
    else:
        return seed + "_" + "_".join(map(str, name_or_coeffs))

def variable(ring: Ring, name_or_coeffs: VariableName, seed: str = "v") -> Variable:
    """ Get variable of a ring from it's name or sequence of coefficients (typically a weight) """
    return variables(ring, (name_or_coeffs,), seed)[0]
    
def variables(ring: Ring, name_or_coeffs: Iterable[VariableName], seed: str = "v") -> tuple[Variable, ...]:
    """ Get multiple variables of a ring from it's name or sequence of coefficients (typically a weight) """
    ring_gens = ring.gens_dict()
    return tuple(ring_gens[variable_name(nc, seed)] for nc in name_or_coeffs )
