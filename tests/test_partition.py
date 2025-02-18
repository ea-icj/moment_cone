import unittest

from cone.partition import Partition

class TestPartition(unittest.TestCase):

    def test_interface(self):
        pe = Partition(())
        self.assertEqual(len(pe), 0)
        self.assertEqual(pe[3], 0)
        self.assertEqual(list(pe), [])

        p = Partition([4, 4, 3, 2, 1, 0, 0])
        self.assertEqual(len(p), 5)

        self.assertEqual(p[0], 4)
        self.assertEqual(p[3], 2)
        self.assertEqual(p[10], 0)

        self.assertEqual(list(p), [4, 4, 3, 2, 1])

        self.assertEqual(p.pad(6), (4, 4, 3, 2, 1, 0))
        self.assertEqual(p.pad(5), (4, 4, 3, 2, 1))
        with self.assertRaises(Exception):
            p.pad(4)

        self.assertTrue(Partition((3, 2, 1)) <= Partition((2, 2, 1)))
        self.assertTrue(Partition((3, 2, 1)) < Partition((2, 2, 1)))
        self.assertFalse(Partition((3, 2, 1)) >= Partition((2, 2, 1)))
        self.assertFalse(Partition((3, 2, 1)) > Partition((2, 2, 1)))

    def test_non_decreasing(self):
        with self.assertRaises(Exception):
            p = Partition([4, 4, 3, 2, 3, 0])

    def test_negative(self):
        with self.assertRaises(Exception):
            p = Partition([4, 4, 3, 2, -1])

    def test_all_of_integer(self):
        partitions = Partition.all_for_integer(4)
        ref = ((4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1))
        for p, r in zip(partitions, ref):
            self.assertEqual(p, Partition(r))

        partitions_len3_ref = [p for p in Partition.all_for_integer(10) if len(p) <= 3]
        partitions_len3 = list(Partition.all_for_integer(10, max_length=3))
        self.assertEqual(partitions_len3, partitions_len3_ref)

    def test_all_of_length(self):
        partitions = Partition.all_of_length(2, 3)
        ref = ((3, 3), (3, 2), (3, 1), (3,), (2, 2), (2, 1), (2,), (1, 1), (1,), ())
        for p, r in zip(partitions, ref):
            self.assertEqual(p, Partition(r))

