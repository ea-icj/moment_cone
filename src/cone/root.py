from .typing import *
from .dimension import Dimension
from dataclasses import dataclass

@dataclass(frozen=True, slots=True)
class Root:
    """ Root element for tau """
    k: int
    i: int
    j: int

    @property
    def is_in_U(self) -> bool:
        """ Check if this root is in U """
        return self.i < self.j


if False: # TODO
    def all_weights_U(d: Dimension) -> Iterable[WeightU]:
        """ Returns all possible weights of U for a given sequence of length """
        # TODO: verify that the actual definition of this weight are so that i < j
        for k, l in enumerate(d):
            for i, j in itertools.combinations(range(l), 2):
                yield k, i, j
