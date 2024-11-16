import unittest

from cone.weight import Weight

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
        d = (4, 2, 3)

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

        

        