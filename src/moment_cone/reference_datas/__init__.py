"""
Load all reference inequalities from reference articles
or from reference executions of this library.
"""
__all__ = (
    'get_reference_ineqs',
    'compare',
    'compare_to_reference',
    'compare_ineq_mod_sym_dim',
)

import os
import sys
import importlib.util
import pathlib
import pkgutil

from ..representation import *
from ..inequality import Inequality
from ..typing import *

from .comparisons import compare_to_reference, compare, compare_ineq_mod_sym_dim

RepresentationKeys: TypeAlias = tuple[
    tuple[int, ...], # linear group
    Optional[int],   # particle count
]

reference_ineqs: dict[RepresentationKeys, tuple[str, list[Inequality]]] = {}

def get_representation_keys(V: Representation) -> RepresentationKeys:
    """ Returns essential properties of the representation
    
    that are used to identify reference inequalities.
    """
    if isinstance(V, ParticleRepresentation):
        return tuple(V.G), V.particle_cnt
    else:
        return tuple(V.G), None
    
def get_reference_ineqs(V: Representation) -> tuple[str, list[Inequality]]:
    """ Returns reference inequalities for the given Representation 
    
    Throw KeyError if no reference exists for this representation.
    """
    return reference_ineqs[get_representation_keys(V)]


# Load all reference inequalities found in this folder
for info in pkgutil.iter_modules(__path__):
    if not info.name.startswith("ineq_"):
        continue
    module = importlib.import_module(f".{info.name}", package=__name__)
    V = cast(Representation, module.V)
    reference_ineqs[get_representation_keys(V)] = module.source, module.inequalities
