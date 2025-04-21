import unittest

from moment_cone.kronecker import *
from moment_cone.partition import Partition

class TestKronecker(unittest.TestCase):
    def test_kronecker(self) -> None:
        partitions_as_tuple = (3, 1, 1, 1, 1, 1), (4, 3, 1), (2, 2, 2, 2), (5, 1, 1, 1), (4, 2, 1,1)
        partitions = tuple(Partition(p) for p in partitions_as_tuple)

        a = Partition(3, 2, 1)
        b = Partition(2, 2, 2)
        c = Partition(2, 2, 1, 1)

        for cls in KroneckerCoefficient, KroneckerCoefficientCache, KroneckerCoefficientMLCache:
            kc = cls()
            self.assertEqual(kc(partitions), 1607)
            self.assertEqual(kc(partitions[:3]), 1)
            self.assertEqual(kc(partitions[:2]), 0)
            self.assertEqual(kc(partitions[:1]), 1) # FIXME: verify this
            self.assertEqual(kc((partitions[0], partitions[0])), 1)

            product = kc.product((a, b, c))
            self.assertEqual(len(product), 11)
            self.assertEqual(product[Partition(6)], 1)
            self.assertEqual(product[Partition(5, 1)], 5)
            self.assertEqual(product[Partition(4, 2)], 9)
            self.assertEqual(product[Partition(4, 1, 1)], 10)
            self.assertEqual(product[Partition(3, 3)], 5)
            self.assertEqual(product[Partition(3, 2, 1)], 16)
            self.assertEqual(product[Partition(3, 1, 1, 1)], 10)
            self.assertEqual(product[Partition(2, 2, 2)], 5)
            self.assertEqual(product[Partition(2, 2, 1, 1)], 9)
            self.assertEqual(product[Partition(2, 1, 1, 1, 1)], 5)
            self.assertEqual(product[Partition(1, 1, 1, 1, 1, 1)], 1)

            product = kc.product((Partition(3, 2, 1),))
            self.assertEqual(len(product), 1)
            self.assertEqual(product[Partition(3, 2, 1)], 1)

            product = kc.product(())
            self.assertEqual(product, dict())
