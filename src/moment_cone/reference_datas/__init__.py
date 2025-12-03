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

reference_ineqs: dict[RepresentationKeys, dict[str, list[Inequality]]] = {}

def get_representation_keys(V: Representation) -> RepresentationKeys:
    """ Returns essential properties of the representation
    
    that are used to identify reference inequalities.
    """
    if isinstance(V, ParticleRepresentation):
        return tuple(V.G), V.particle_cnt
    else:
        return tuple(V.G), None
    
def get_reference_ineqs(V: Representation, source: Optional[str] = None) -> tuple[str, list[Inequality]]:
    """ Returns reference inequalities for the given Representation and optional source
    
    Throw KeyError if no reference exists for this representation.
    """
    V_keys = get_representation_keys(V)
    if V_keys not in reference_ineqs:
        load_all_ineq_for_repr(V)
    ineq_by_source = reference_ineqs[V_keys]

    if source is None:
        return next(iter(ineq_by_source.items()))
    else:
        return source, ineq_by_source[source.lower()]

def load_ineq_from_name(name: str) -> None:
    """ Load inequalities for a given submodule name from reference_datas folder """
    module = importlib.import_module(f".{name}", package=__name__)
    V = cast(Representation, module.V)
    source = cast(str, module.source)
    inequalities = cast(list[Inequality], module.inequalities)
    ineq_by_source = reference_ineqs.setdefault(get_representation_keys(V), dict())
    ineq_by_source[source.lower()] = inequalities

def load_all_ineq_for_repr(V: Representation) -> None:
    """ Load inequalities for all sources of a given representation """
    from ..export import generate_file_name
    base_name = generate_file_name(V)[:-1]
    for info in pkgutil.iter_modules(__path__):
        # Not specific enough but OK
        if info.name.lower().startswith(base_name.lower()):
            load_ineq_from_name(info.name)

def load_all_ineq() -> None:
    """ Load all reference inequalities from reference_datas folder """
    for info in pkgutil.iter_modules(__path__):
        if info.name.lower().startswith("ineq_"):
            load_ineq_from_name(info.name)
