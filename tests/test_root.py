import unittest

from cone.root import Root
from cone.dimension import Dimension

class TestRoot(unittest.TestCase):
    def test_interface(self) -> None:
        r1 = Root(1, 2, 3)
        r2 = Root(1, 3, 2)

        self.assertEqual(r1.k, 1)
        self.assertEqual(r1.i, 2)
        self.assertEqual(r1.j, 3)

        self.assertTrue(r1.is_in_U)
        self.assertFalse(r2.is_in_U)

    def test_all(self) -> None:
        d = Dimension((2, 3, 4))
        all_r = tuple(Root.all_of_U(d))

        self.assertEqual(len(all_r), 10)
        self.assertEqual(all_r[0], Root(0, 0, 1))
        self.assertEqual(all_r[3], Root(1, 1, 2))
        self.assertEqual(all_r[5], Root(2, 0, 2))
        self.assertEqual(all_r[9], Root(2, 2, 3))