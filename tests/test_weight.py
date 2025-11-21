import unittest

from moment_cone.weight import WeightAsList
from moment_cone.linear_group import LinearGroup
from moment_cone.representation import KroneckerRepresentation

class TestWeightAsList(unittest.TestCase):

    def test_interface(self) -> None:
        G = LinearGroup((2, 3, 6, 1))
        wt = (1, 2, 4, 0)
        w = WeightAsList(G, wt)

        self.assertEqual(len(w), len(wt))
        self.assertEqual(w.as_list[1], wt[1])
        self.assertEqual(tuple(w), wt)

    def test_comparison(self) -> None:
        G = LinearGroup((2, 4, 6, 2))
        w1 = WeightAsList(G, [1, 2, 4, 0], index=1) # Dummy index
        w2 = WeightAsList(G, [1, 2, 5, 0], index=2) # Dummy index
        w3 = WeightAsList(G, [1, 3, 3, 0], index=3) # Dummy index
        w4 = WeightAsList(G, [1, 2, 4, 0], index=4) # Dummy index

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

    def test_all_index(self) -> None:
        import functools
        import operator
        G = LinearGroup((4, 2, 3, 1))
        V = KroneckerRepresentation(G)

        all_weights = list(V.all_weights)

        self.assertEqual(len(all_weights), functools.reduce(operator.mul, G))
        self.assertEqual(all_weights[8], WeightAsList(G, [1, 0, 2, 0]))
        self.assertEqual(all_weights[-1], WeightAsList(G, [di - 1 for di in G]))

        for i, w in enumerate(all_weights):
            self.assertEqual(i, w.index)
            # self.assertEqual(w, Weight.from_index(d, i)) # FIXME
            self.assertEqual(i, V.index_of_weight(w))
            self.assertEqual(i, V.index_of_weight(w, use_internal_index=False))
            self.assertEqual(i, V.index_of_weight(WeightAsList(G, list(w))))

    
    def test_all_mod_sym_dim(self) -> None:
        from itertools import pairwise

        G = LinearGroup((4, 4, 4, 2, 2))
        V = KroneckerRepresentation(G)

        # FIXME
        """
        mod_weights = list(Weight.all_mod_sym_dim(d))

        mod_weights_ref = []
        for w in Weight.all(d):
            wl = tuple(w)
            if all(a >= b for a, b in pairwise(wl[:3])) and all(a >= b for a, b in pairwise(wl[3:])):
                mod_weights_ref.append(w)

        self.assertEqual(len(mod_weights), len(mod_weights_ref))
        self.assertTrue(all(w in mod_weights_ref for w in mod_weights))

        self.assertEqual(mod_weights[0], WeightAsList((3, 3, 3, 1, 1)))
        self.assertEqual(mod_weights[4], WeightAsList((3, 3, 2, 1, 0)))
        self.assertEqual(mod_weights[-1], WeightAsList((0, 0, 0, 0, 0)))
        """

    def test_orbit(self) -> None:
        from moment_cone.utils import symmetries
        G = LinearGroup((5, 5, 5, 3, 3, 3, 5))
        p = WeightAsList(G, (2, 2, 4, 1, 2, 1, 4))
        orbits = list(p.orbit_symmetries(symmetries(G)))
        self.assertEqual(len(orbits), 9)
        self.assertEqual(orbits[0], WeightAsList(G, (2, 2, 4, 1, 1, 2, 4)))
        self.assertEqual(orbits[1], WeightAsList(G, (2, 2, 4, 1, 2, 1, 4)))
        self.assertEqual(orbits[2], WeightAsList(G, (2, 2, 4, 2, 1, 1, 4)))
        self.assertEqual(orbits[3], WeightAsList(G, (2, 4, 2, 1, 1, 2, 4)))
        self.assertEqual(orbits[4], WeightAsList(G, (2, 4, 2, 1, 2, 1, 4)))
        self.assertEqual(orbits[5], WeightAsList(G, (2, 4, 2, 2, 1, 1, 4)))
        self.assertEqual(orbits[6], WeightAsList(G, (4, 2, 2, 1, 1, 2, 4)))
        self.assertEqual(orbits[7], WeightAsList(G, (4, 2, 2, 1, 2, 1, 4)))
        self.assertEqual(orbits[8], WeightAsList(G, (4, 2, 2, 2, 1, 1, 4)))
        