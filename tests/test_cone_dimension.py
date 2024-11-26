import unittest
import random

from sage.all import ZZ, QQ, I # type: ignore

from cone.dimension import Dimension
from cone.weight import Weight
import cone.cone_dimension as cd

class TestPermutation(unittest.TestCase):
    def test_pointv(self) -> None:
        d = Dimension((4, 3, 2))
        N = 10
        pds = random.sample(tuple(Weight.all(d)), N)

        for ring in ZZ, QQ, QQ[I]:
            v = cd.PointV(pds, d, ring, upper=10)
            self.assertTrue(all(
                v[i] == 0 for i in range(d.prod)
                if Weight.from_index(d, i) not in pds
            ))
            self.assertTrue(all(
                1 <= v[chi.index] <= 10
                for chi in pds
            ))