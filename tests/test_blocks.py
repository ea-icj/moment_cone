import unittest

from cone.blocks import Blocks

class TestBlocks(unittest.TestCase):

    def test_mutable(self) -> None:
        b = Blocks.from_flatten([1, 2, 3, 4], (2, 1, 1))

        self.assertEqual(repr(b), "1 2 | 3 | 4")
        self.assertFalse(b.is_frozen)
        self.assertEqual(b.sizes, (2, 1, 1))
        self.assertEqual(len(b), 3)
        self.assertEqual(tuple(b.blocks), ([1, 2], [3], [4]))
        self.assertEqual(tuple(b), ([1, 2], [3], [4]))
        self.assertEqual(b.flatten, [1, 2, 3, 4])

        self.assertEqual(b[0], [1, 2])
        self.assertEqual(b[1], [3])
        self.assertEqual(b[0, 1], 2)
        self.assertEqual(tuple(b[1:]), ([3], [4]))

        b[0, 1] = 8
        b[1] = [3, 4, 5]

        self.assertEqual(repr(b), "1 8 | 3 4 5 | 4")
        self.assertFalse(b.is_frozen)
        self.assertEqual(b.sizes, (2, 3, 1))
        self.assertEqual(len(b), 3)
        self.assertEqual(tuple(b.blocks), ([1, 8], [3, 4, 5], [4]))
        self.assertEqual(tuple(b), ([1, 8], [3, 4, 5], [4]))
        self.assertEqual(b.flatten, [1, 8, 3, 4, 5, 4])

        self.assertEqual(b[0], [1, 8])
        self.assertEqual(b[1], [3, 4, 5])
        self.assertEqual(b[0, 1], 8)
        self.assertEqual(tuple(b[1:]), ([3, 4, 5], [4]))

        b2 = Blocks.from_flatten([1, 8, 3, 4, 5, 4], (2, 3, 1))
        self.assertEqual(b, b2)
        self.assertNotEqual(hash(b), hash(b2))

    def test_immutable(self) -> None:
        b = Blocks.from_flatten((1, 2, 3, 4), (2, 1, 1))

        self.assertEqual(repr(b), "1 2 | 3 | 4")
        self.assertTrue(b.is_frozen)
        self.assertEqual(b.sizes, (2, 1, 1))
        self.assertEqual(len(b), 3)
        self.assertEqual(tuple(b.blocks), ((1, 2), (3,), (4,)))
        self.assertEqual(tuple(b), ((1, 2), (3,), (4,)))
        self.assertEqual(b.flatten, (1, 2, 3, 4))

        self.assertEqual(b[0], (1, 2))
        self.assertEqual(b[1], (3,))
        self.assertEqual(b[0, 1], 2)
        self.assertEqual(tuple(b[1:]), ((3,), (4,)))

        with self.assertRaises(Exception):
            b[0, 1] = 8
        
        with self.assertRaises(Exception):
            b[1] = [3, 4, 5]

        b2 = Blocks.from_flatten((1, 2, 3, 4), (2, 1, 1))
        self.assertEqual(b, b2)
        self.assertEqual(hash(b), hash(b2))

        b3 = Blocks.from_blocks(((1, 2), (3,), (4,)))
        self.assertEqual(b, b3)
        self.assertEqual(hash(b), hash(b3))

        b4 = Blocks.from_flatten([1, 2, 3, 4], (2, 1, 1)).freeze()
        self.assertEqual(b, b4)
        self.assertEqual(hash(b), hash(b4))

    def test_sage_immutable(self) -> None:
        from sage.all import Integer # type: ignore
        from typing import Iterable
        def to_sage(*s: int) -> tuple[Integer]:
            return tuple(map(Integer, s))
        
        b = Blocks.from_flatten(to_sage(1, 2, 3, 4), to_sage(2, 1, 1))

        self.assertEqual(repr(b), "1 2 | 3 | 4")
        self.assertTrue(b.is_frozen)
        self.assertEqual(b.sizes, (2, 1, 1))
        self.assertEqual(len(b), 3)
        self.assertEqual(tuple(b.blocks), ((1, 2), (3,), (4,)))
        self.assertEqual(tuple(b), ((1, 2), (3,), (4,)))
        self.assertEqual(b.flatten, (1, 2, 3, 4))

        self.assertEqual(b[Integer(0)], (1, 2))
        self.assertEqual(b[Integer(1)], (3,))
        self.assertEqual(b[Integer(0), Integer(1)], 2)
        self.assertEqual(tuple(b[Integer(1):]), ((3,), (4,)))

        with self.assertRaises(Exception):
            b[Integer(0), Integer(1)] = 8
        
        with self.assertRaises(Exception):
            b[Integer(1)] = [3, 4, 5]

        b2 = Blocks.from_flatten(to_sage(1, 2, 3, 4), to_sage(2, 1, 1))
        self.assertEqual(b, b2)
        self.assertEqual(hash(b), hash(b2))

        b3 = Blocks.from_blocks((to_sage(1, 2), to_sage(3), to_sage(4)))
        self.assertEqual(b, b3)
        self.assertEqual(hash(b), hash(b3))

        b4 = Blocks.from_flatten(list(to_sage(1, 2, 3, 4)), to_sage(2, 1, 1)).freeze()
        self.assertEqual(b, b4)
        self.assertEqual(hash(b), hash(b4))
