import unittest

from cone.weight import Weight
from cone.dimension import Dimension

class TestWeight(unittest.TestCase):

    def test_interface(self):
        wt = (1, 2, 4, 0)
        w = Weight(wt)

        self.assertEqual(len(w), len(wt))
        self.assertEqual(w[1], wt[1])
        self.assertEqual(tuple(w), wt)

    def test_comparison(self):
        w1 = Weight([1, 2, 4, 0], 1) # Dummy index
        w2 = Weight([1, 2, 5, 0], 2) # Dummy index
        w3 = Weight([1, 3, 3, 0], 3) # Dummy index
        w4 = Weight([1, 2, 4, 0], 4) # Dummy index

        self.assertTrue(w1 == w4) # Equal even if index is different

        self.assertTrue(w1 >= w2)
        self.assertFalse(w1 <= w2)
        self.assertTrue(w2 <= w1)
        self.assertFalse(w2 >= w1)

        self.assertTrue(w1 >= w2)
        self.assertFalse(w1 <= w2)
        self.assertTrue(w2 <= w1)
        self.assertFalse(w2 >= w1)

        # Partial ordering
        self.assertFalse(w1 <= w3)
        self.assertFalse(w1 >= w3)
        self.assertFalse(w3 >= w1)
        self.assertFalse(w3 <= w1)

    def test_all_index(self):
        import functools
        import operator
        d = Dimension((4, 2, 3))

        all_weights = list(Weight.all(d))

        self.assertEqual(len(all_weights), functools.reduce(operator.mul, d))
        self.assertEqual(all_weights[8], Weight([1, 0, 2]))
        self.assertEqual(all_weights[-1], Weight([di - 1 for di in d]))

        for i, w in enumerate(all_weights):
            self.assertEqual(i, w.index)
            self.assertEqual(w, Weight.from_index(d, i))
            self.assertEqual(i, w.index_in(d))
            self.assertEqual(i, w.index_in(d, use_internal_index=False))
            self.assertEqual(i, Weight(list(w)).index_in(d))

    
    def test_all_mod_sym_dim(self):
        from itertools import pairwise

        d = Dimension((4, 4, 4, 2, 2))
        mod_weights = list(Weight.all_mod_sym_dim(d))

        mod_weights_ref = []
        for w in Weight.all(d):
            wl = tuple(w)
            if all(a >= b for a, b in pairwise(wl[:3])) and all(a >= b for a, b in pairwise(wl[3:])):
                mod_weights_ref.append(w)

        self.assertEqual(len(mod_weights), len(mod_weights_ref))
        self.assertTrue(all(w in mod_weights_ref for w in mod_weights))

        self.assertEqual(mod_weights[0], Weight((3, 3, 3, 1, 1)))
        self.assertEqual(mod_weights[4], Weight((3, 3, 2, 1, 0)))
        self.assertEqual(mod_weights[-1], Weight((0, 0, 0, 0, 0)))


    def test_orbit(self) -> None:
        p = Weight((2, 2, 4, 1, 2, 1, 4))
        orbits = list(p.orbit_symmetries((3, 3, 1)))
        self.assertEqual(len(orbits), 9)
        self.assertEqual(orbits[0], Weight((2, 2, 4, 1, 1, 2, 4)))
        self.assertEqual(orbits[1], Weight((2, 2, 4, 1, 2, 1, 4)))
        self.assertEqual(orbits[2], Weight((2, 2, 4, 2, 1, 1, 4)))
        self.assertEqual(orbits[3], Weight((2, 4, 2, 1, 1, 2, 4)))
        self.assertEqual(orbits[4], Weight((2, 4, 2, 1, 2, 1, 4)))
        self.assertEqual(orbits[5], Weight((2, 4, 2, 2, 1, 1, 4)))
        self.assertEqual(orbits[6], Weight((4, 2, 2, 1, 1, 2, 4)))
        self.assertEqual(orbits[7], Weight((4, 2, 2, 1, 2, 1, 4)))
        self.assertEqual(orbits[8], Weight((4, 2, 2, 2, 1, 1, 4)))
        