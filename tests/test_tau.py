import unittest

from cone.tau import Tau, ReducedTau
from cone.blocks import Blocks
from cone.dimension import Dimension
from cone.root import Root

class TestTau(unittest.TestCase):

    def test_init(self):
        d = Dimension((3, 2, 4))
        tau_columns = Tau(((3, 2, 1), (2, 3), (1, 4, 5, 3)), 3)
        tau_flatten = Tau.from_flatten((3, 3, 2, 1, 2, 3, 1, 4, 5, 3), d)
        
        blocks = Blocks.from_blocks(((3, 2, 1), (2, 3), (1, 4, 5, 3)))
        tau_matrix = Tau(blocks, 3)

        self.assertTrue(tau_columns.ccomponent == tau_flatten.ccomponent == tau_matrix.ccomponent)
        self.assertTrue(tau_columns.components == tau_flatten.components == tau_matrix.components)

    def test_base_interface(self):
        tau = Tau(((3, 2, 1), (2, 3), (1, 4, 5, 3)), 3)
        self.assertEqual(len(tau), 3)
        self.assertEqual(tau.d, (3, 2, 4))
        self.assertEqual(tau.components, ((3, 2, 1), (2, 3), (1, 4, 5, 3)))
        self.assertEqual(tuple(tau.flattened), (3, 3, 2, 1, 2, 3, 1, 4, 5, 3))
        self.assertEqual(repr(tau), "3 | 3 2 1 | 2 3 | 1 4 5 3")

    def test_dot(self):
        # TODO
        pass

    def test_regular(self):
        tau1 = Tau(((3, 2, 2), (4, 2, 1), (3, 2)))
        tau2 = Tau(((3, 2, 2), (4, 2, 3), (3, 2)))
        tau3 = Tau(((3, 2, 1), (4, 2, 1), (3, 2)))

        self.assertTrue(tau1.is_dominant)
        self.assertFalse(tau2.is_dominant)
        self.assertTrue(tau3.is_dominant)

        self.assertFalse(tau1.is_dom_reg)
        self.assertTrue(tau3.is_dom_reg)
        
    def test_representative(self):
        d = Dimension((3, 3, 3))
        tau = Tau.from_flatten([1, 3, 2, 1, 4, 1, 2, 5, 3, 1], d)
        self.assertEqual(repr(tau.sl_representative), "25 | 3 0 -3 | 5 -4 -1 | 6 0 -6")
        self.assertEqual(repr(tau.end0_representative), "5 | 2 1 0 | 2 -1 0 | 4 2 0")

    def test_pos_weights(self):
        # TODO
        pass

    def test_pos_roots(self):
        tau = Tau(((6, 2), (1, 4, 1), (2, 5, 3, 1)), 1)
        pr = tau.positive_roots
        self.assertEqual(sorted(pr.keys()), [1, 2, 3, 4])
        self.assertEqual(pr[1], [Root(2, 0, 3)])
        self.assertEqual(pr[2], [Root(2, 1, 2), Root(2, 2, 3)])
        self.assertEqual(pr[3], [Root(1, 1, 2)])
        self.assertEqual(pr[4], [Root(0, 0, 1), Root(2, 1, 3)])

    def test_sort_mod_sd(self):
        d = Dimension((2, 2, 2, 1, 1, 1))
        tau = Tau.from_flatten([1,6,2,1,4,1,2,5,3,1], d)
        self.assertEqual(repr(tau.sort_mod_sym_dim), "1 | 1 2 | 1 4 | 6 2 | 1 | 3 | 5")

    def test_reduced_tau(self):
        tau = Tau(((2, 2, 3), (4, 3, 3, 2, 2, 2, 1), (5, 2, 2)), 3)
        red_tau = tau.reduced

        self.assertEqual(len(red_tau), 3)
        self.assertEqual(red_tau.small_d, (2, 4, 2))
        self.assertEqual(red_tau.ccomponent, 3)
        self.assertEqual(tuple(red_tau.mult.blocks), ((2, 1), (1, 2, 3, 1), (1, 2)))
        self.assertEqual(tuple(red_tau.values.blocks), ((2, 3), (4, 3, 2, 1), (5, 2)))

        self.assertEqual(red_tau[1, 2], (2, 2))
        self.assertEqual(red_tau[0, 1], (4, 1))

    def test_pzero(self):
        # TODO
        pass

    def test_hash(self):
        d = Dimension((2, 3, 4))
        tau1 = Tau.from_flatten([6,2,1,4,1,2,5,3,1], d)
        tau2 = Tau.from_flatten([6,2,1,4,1,2,5,3,1], d)
        self.assertEqual(tau1, tau2)
        self.assertEqual(hash(tau1), hash(tau2))

        rtau1 = tau1.reduced
        rtau2 = tau2.reduced
        self.assertEqual(rtau1, rtau2)
        self.assertEqual(hash(rtau1), hash(rtau2))
