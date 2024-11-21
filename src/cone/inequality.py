from .typing import *
from .tau import Tau
from .permutation import Permutation

class Inequality:
    """ An inequality composed of a Tau and a permutation """
    tau: Tau
    w: tuple[Permutation, ...]
    wtau: Tau # each w_k applied to each column of tau

    def __init__(self, tau: Tau, w: Iterable[Permutation]):
        self.tau = tau
        self.w = tuple(w)
        assert len(tau.d) == len(self.w)
        self.wtau = Tau(tuple(wk(ck) for wk, ck in zip(self.w, tau.components)), tau.ccomponent)

    def __repr__(self) -> str:
        return \
            f"Inequality(tau  = {self.tau},\n" + \
             "           w    =     " + " | ".join(" ".join(map(str, wk)) for wk in self.w) + ",\n" + \
            f"           wtau = {self.wtau})"
