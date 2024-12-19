import unittest
import random

from sage.all import ZZ, QQ, I # type: ignore

from cone.dimension import Dimension
from cone.weight import Weight
import cone.vector_chooser as cd

class TestPermutation(unittest.TestCase):
    def test_pointv(self) -> None:
        d = Dimension((4, 3, 2))
        N = 10
        pds = random.sample(tuple(Weight.all(d)), N)

        for ring in d.Q, d.QI, d.QV, d.QIV, d.QV2:
            v = cd.point_vect(pds, d, ring, bounds=(4, 10))

            self.assertTrue(all(
                v[i] == 0 for i in range(d.dimV)
                if Weight.from_index(d, i) not in pds
            ))

            if ring is d.Q:
                self.assertTrue(all(
                    4 <= v[chi.index] <= 10
                    for chi in pds
                ))
            elif ring is d.QI:
                self.assertTrue(all(
                    4 <= v[chi.index].real() <= 10
                    and
                    4 <= v[chi.index].imag() <= 10
                    for chi in pds
                ))
            else:
                pass # TODO
