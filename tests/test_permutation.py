import unittest

from cone.permutation import Permutation
from cone.dimension import Dimension

class TestPermutation(unittest.TestCase):
    def test_interface(self) -> None:
        # (1,20,16,12,8,4,21,17,13,9,5)(2,18,14,10,6)(3,19,15,11,7)
        # with cycle length 11, 5, 5
        wp = (19, 17, 18, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
        p = Permutation(wp)
        self.assertEqual(len(p), 21)
        self.assertEqual(p.n, 21)
        self.assertEqual(p.length, 70)
        self.assertEqual(p(range(21)), wp)
        self.assertEqual(p(p(range(21))), (15, 13, 14, 16, 19, 17, 18, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

        s = tuple(range(21))
        for _ in range(55): # lcm(11, 5, 5) = 55
            s = p(s)
        self.assertEqual(s, tuple(range(21)))

        self.assertEqual(p.inverse(p(range(21))), tuple(range(21)))
        self.assertEqual(p, p.inverse.inverse)

    def test_all(self) -> None:
        pl = list(Permutation.all(3))
        self.assertEqual(len(pl), 6)
        self.assertEqual(pl[0], Permutation((0, 1, 2)))
        self.assertEqual(pl[1], Permutation((0, 2, 1)))
        self.assertEqual(pl[2], Permutation((1, 0, 2)))
        self.assertEqual(pl[3], Permutation((1, 2, 0)))
        self.assertEqual(pl[4], Permutation((2, 0, 1)))
        self.assertEqual(pl[5], Permutation((2, 1, 0)))

    def test_all_of_length(self) -> None:
        pl = list(Permutation.all_of_length(4, 3))
        self.assertEqual(len(pl), 6)
        self.assertEqual(pl[0], Permutation((0, 3, 2, 1)))
        self.assertEqual(pl[1], Permutation((1, 2, 3, 0)))
        self.assertEqual(pl[2], Permutation((1, 3, 0, 2)))
        self.assertEqual(pl[3], Permutation((2, 0, 3, 1)))
        self.assertEqual(pl[4], Permutation((2, 1, 0, 3)))
        self.assertEqual(pl[5], Permutation((3, 0, 1, 2)))

    def test_others(self) -> None:
        p = Permutation((3, 2, 1, 4, 3, 5, 4, 3, 2))
        
        d1 = Dimension((3, 2, 4))
        self.assertTrue(p.is_min_rep(d1))

        d2 = Dimension((2, 3, 4))
        self.assertFalse(p.is_min_rep(d2))

    def test_permutationby_length(self) -> None:
        # TODO
        pass