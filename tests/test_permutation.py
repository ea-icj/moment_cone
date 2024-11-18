import unittest

from cone.permutation import Permutation

class TestPermutation(unittest.TestCase):
    def test_interface(self) -> None:
        # (1,20,16,12,8,4,21,17,13,9,5)(2,18,14,10,6)(3,19,15,11,7)
        # with cycle length 11, 5, 5
        wp = (19, 17, 18, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
        p = Permutation(wp)
        self.assertEqual(len(p), 21)
        self.assertEqual(p.length, 70)
        self.assertEqual(p(range(21)), wp)
        self.assertEqual(p(p(range(21))), (15, 13, 14, 16, 19, 17, 18, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

        s = tuple(range(21))
        for _ in range(55): # lcm(11, 5, 5) = 55
            s = p(s)
        self.assertEqual(s, tuple(range(21)))

        self.assertEqual(p.inverse(p(range(21))), tuple(range(21)))

    def test_all(self) -> None:
        # TODO: testing generators
        pass
