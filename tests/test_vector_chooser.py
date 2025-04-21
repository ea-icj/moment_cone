import unittest
import random

from sage.all import ZZ, QQ, I # type: ignore

from moment_cone.linear_group import LinearGroup
from moment_cone.weight import Weight
from moment_cone.representation import KroneckerRepresentation
import moment_cone.vector_chooser as cd

class TestPermutation(unittest.TestCase):
    def test_pointv(self) -> None:
        G = LinearGroup((4, 3, 2))
        V = KroneckerRepresentation(G)
        N = 10
        pds = random.sample(tuple(V.all_weights), N)

        for ring in V.Q, V.QI, V.QV, V.QIV:
            v = cd.point_vect(pds, V, ring, bounds=(4, 10))

            # FIXME
            """
            self.assertTrue(all(
                v[i] == 0 for i in range(d.dimV)
                if Weight.from_index(d, i) not in pds
            ))
            """

            if ring is V.Q:
                self.assertTrue(all(
                    4 <= v[V.index_of_weight(chi)] <= 10
                    for chi in pds
                ))
            elif ring is V.QI:
                self.assertTrue(all(
                    4 <= v[V.index_of_weight(chi)].real() <= 10
                    and
                    4 <= v[V.index_of_weight(chi)].imag() <= 10
                    for chi in pds
                ))
            else:
                pass # TODO
