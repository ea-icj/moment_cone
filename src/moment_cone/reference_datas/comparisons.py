__all__ = (
    "Comparison",
    "compare",
    "compare_ineq_mod_sym_dim",
    "compare_to_reference",
)

from functools import cached_property

from ..tau import *
from ..representation import *
from ..inequality import *
from ..typing import *

class Comparison(Generic[T]):
    """ Compare two lists of elements """
    set1: set[T]
    set2: set[T]
    element_name: str
    source1: str
    source2: str

    def __init__(
        self,
        list1: Iterable[T], #: first list
        list2: Iterable[T], #: second list
        element_name: Optional[str] = None, #: display name of the elements
        source1: str = "1", #: origin of the first list
        source2: str = "2", #: origin of the second list
    ):
        self.set1 = set(list1)
        self.set2 = set(list2)

        if element_name is None:
            from itertools import chain
            try:
                an_element = next(chain(self.set1, self.set2))
            except StopIteration:
                element_name = "element"
            else:
                element_name = type(an_element).__name__
        self.element_name = element_name

        self.source1 = source1
        self.source2 = source2

    def __bool__(self) -> bool:
        """ Returns true if the two lists are identical """
        return self.set1 == self.set2
    
    @cached_property
    def intersection(self) -> set[T]:
        return self.set1 & self.set2
    
    @cached_property
    def only1(self) -> set[T]:
        return self.set1 - self.set2

    @cached_property
    def only2(self) -> set[T]:
        return self.set2 - self.set1

    def __repr__(self) -> str:
        from textwrap import dedent
        return dedent(f"""
        The two lists of {self.element_name} are {"identical" if self else "different"}:
        {len(self.set1)} {self.element_name} from first list of {self.source1} vs {len(self.set2)} {self.element_name} from second list of {self.source2}
        {len(self.intersection)} {self.element_name} are in both lists
        {len(self.only1)} {self.element_name} appearing only in the first list of {self.source1}
        {len(self.only2)} {self.element_name} appearing only in the second list of {self.source2}
        """).strip()


#TODO check unitary tests
def compare(
        list1: Iterable[T],
        list2: Iterable[T],
        element_name: Optional[str] = None,
        source1: str = "1",
        source2: str = "2",
    ) -> Comparison[T]:
    """ Compares two lists. returns elements from the first one that don't belong to the second one and conversely.
    elements have to be hashable (e.g. tuples) so that set(list1) works
    optional arguments comments are here to caracterie the nature of elements (comment0) and particularities of the compared lists (cf compareK_ineq and compareVW_ineq below)

    >>> compare([(1,3,6),(2,4,9)],[(1,3,6),(12,43,24),(1,3,6)])
    2 elements from list 1  vs  2 elements from list 2
    1 elements in both lists
    1 elements appearing only in list 1 {(2, 4, 9)}
    1 elements appearing only in list 2 {(12, 43, 24)}
    [{(2, 4, 9)}, {(12, 43, 24)}]
    """
    return Comparison(list1, list2, element_name, source1, source2)


def compare_ineq_mod_sym_dim(
        list1_ineq: Iterable[Inequality],
        list2_ineq: Iterable[Inequality],
        source1: str = "1",
        source2: str = "2"
    ) -> Comparison[Inequality]:
    list1_ineq = set(list1_ineq)
    list2_ineq = set(list2_ineq)

    from collections import defaultdict
    from itertools import chain
    tau_ineq: dict[Tau, set[Inequality]] = defaultdict(set)
    for ineq in chain(list1_ineq, list2_ineq):
        tau = ineq.wtau.end0_representative.sort_mod_sym_dim
        tau_ineq[tau].add(ineq)

    unique_ineq1: set[Inequality] = set()
    unique_ineq2: set[Inequality] = set()
    for ineqs in tau_ineq.values():
        if ineqs & list1_ineq: unique_ineq1.add(next(iter(ineqs)))
        if ineqs & list2_ineq: unique_ineq2.add(next(iter(ineqs)))

    return compare(unique_ineq1, unique_ineq2, "inequalities (up to S3-sym)", source1, source2)
        

def compare_to_reference(list_ineq: Sequence[Inequality], V: Representation, source: str = "user") -> Comparison[Inequality]:
    """
    list_ineq is a list of Inequalities computed for a certain representation V.
    If exists, it will be compared to a reference list of inequalities (currently only the cases of Klyachko.py for fermions and Vergne_Walter.py for kronecker)
    Other references can be added later
    """
    from . import get_reference_ineqs
    source_ref, reference = get_reference_ineqs(V)
    element_name = "inequalities"
    if isinstance(V, KroneckerRepresentation):
        return compare_ineq_mod_sym_dim(list_ineq, reference, source, source_ref)
    else:
        return compare(list_ineq, reference, element_name, source, source_ref)





