import unittest

from moment_cone.utils import *

class TestUtils(unittest.TestCase):

    def test_monotonicity(self) -> None:
        a = (1, 2, 3, 3, 4)
        b = (4, 4, 3, 1, 1, 0)
        c = (1, 2, 3, 2, 4)
        d = (1, 1, 1, 1, 1)

        self.assertFalse(is_decreasing(a))
        self.assertTrue(is_decreasing(b))
        self.assertFalse(is_decreasing(c))
        self.assertTrue(is_decreasing(d))
    
        self.assertTrue(is_increasing(a))
        self.assertFalse(is_increasing(b))
        self.assertFalse(is_increasing(c))
        self.assertTrue(is_increasing(d))
    
    def test_trim(self) -> None:
        self.assertEqual(trim_zeros((0, 1, 2, 0, 0)), (0, 1, 2))
        self.assertEqual(trim_zeros((0, 1, 2, 0)), (0, 1, 2))
        self.assertEqual(trim_zeros((0, 1, 2)), (0, 1, 2))
        self.assertEqual(trim_zeros((0, 0, 0, 0)), ())

    def test_compression(self) -> None:
        s = (2, 2, 2, 3, 3, 1, 1, 1, 1, 5, 6, 6)
        
        sc = tuple(group_by_block(s))
        self.assertEqual(sc, ((2, 3), (3, 2), (1, 4), (5, 1), (6, 2)))

        sd = tuple(expand_blocks(*zip(*sc)))
        self.assertEqual(sd, s)

        self.assertEqual(tuple(group_by_block(())), ())
        self.assertEqual(tuple(expand_blocks((), ())), ())