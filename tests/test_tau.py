import unittest

from moment_cone.tau import Tau, ReducedTau
from moment_cone.blocks import Blocks
from moment_cone.linear_group import LinearGroup
from moment_cone.root import Root

class TestTau(unittest.TestCase):

    def test_init(self) -> None:
        G = LinearGroup((3, 2, 4, 1))
        tau_columns = Tau(((3, 2, 1), (2, 3), (1, 4, 5, 3), (3,)))
        tau_flatten = Tau.from_flatten((3, 2, 1, 2, 3, 1, 4, 5, 3, 3), G)
        
        blocks = Blocks.from_blocks(((3, 2, 1), (2, 3), (1, 4, 5, 3), (3,)))
        tau_matrix = Tau(blocks)

        self.assertTrue(tau_columns.components == tau_flatten.components == tau_matrix.components)

    def test_base_interface(self) -> None:
        tau = Tau(((3, 2, 1), (2, 3), (1, 4, 5, 3), (3,)))
        self.assertEqual(len(tau), 4)
        self.assertEqual(tau.G, (3, 2, 4, 1))
        self.assertEqual(tau.components, ((3, 2, 1), (2, 3), (1, 4, 5, 3), (3,)))
        self.assertEqual(tuple(tau.flattened), (3, 2, 1, 2, 3, 1, 4, 5, 3, 3))
        self.assertEqual(repr(tau), "3 2 1 | 2 3 | 1 4 5 3 | 3")

    def test_dot(self) -> None:
        # TODO
        pass

    def test_regular(self) -> None:
        tau1 = Tau(((3, 2, 2), (4, 2, 1), (3, 2)))
        tau2 = Tau(((3, 2, 2), (4, 2, 3), (3, 2)))
        tau3 = Tau(((3, 2, 1), (4, 2, 1), (3, 2)))

        self.assertTrue(tau1.is_dominant)
        self.assertFalse(tau2.is_dominant)
        self.assertTrue(tau3.is_dominant)

        self.assertFalse(tau1.is_dom_reg)
        self.assertTrue(tau3.is_dom_reg)
        
    def test_representative(self) -> None:
        G = LinearGroup((3, 3, 3, 1))
        tau = Tau.from_flatten([3, 2, 1, 4, 1, 2, 5, 3, 1, 1], G)
        self.assertEqual(repr(tau.sl_representative), "3 0 -3 | 5 -4 -1 | 6 0 -6 | 25")
        self.assertEqual(repr(tau.end0_representative), "2 1 0 | 2 -1 0 | 4 2 0 | 5")

    def test_pos_weights(self) -> None:
        # TODO
        pass

    def test_pos_grading_rootsU(self) -> None:
        tau = Tau(((6, 2), (1, 4, 1), (2, 5, 3, 1), (1,)))
        pr = tau.grading_rootsU
        self.assertEqual(sorted(pr.keys()), [-3, -1, 0, 1, 2, 3, 4])
        self.assertEqual(pr[-3], [Root(1, 0, 1), Root(2, 0, 1)])
        self.assertEqual(pr[-1], [Root(2, 0, 2)])
        self.assertEqual(pr[0], [Root(1, 0, 2)])
        self.assertEqual(pr[1], [Root(2, 0, 3)])
        self.assertEqual(pr[2], [Root(2, 1, 2), Root(2, 2, 3)])
        self.assertEqual(pr[3], [Root(1, 1, 2)])
        self.assertEqual(pr[4], [Root(0, 0, 1), Root(2, 1, 3)])

    def test_pos_grading_rootsB(self) -> None:
        pass # TODO

    def test_sort_mod_sd(self) -> None:
        G = LinearGroup((1, 2, 2, 2, 1, 1, 1))
        tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 2, 5, 3, 1], G)
        self.assertEqual(repr(tau.sort_mod_sym_dim), "1 | 1 2 | 1 4 | 6 2 | 1 | 3 | 5")

    def test_reduced_tau(self) -> None:
        tau = Tau(((2, 2, 3), (4, 3, 3, 2, 2, 2, 1), (5, 2, 2), (3,)))
        red_tau = tau.reduced

        self.assertEqual(len(red_tau), 4)
        self.assertEqual(red_tau.small_d, (2, 4, 2, 1))
        self.assertEqual(tuple(red_tau.mult.blocks), ((2, 1), (1, 2, 3, 1), (1, 2), (1,)))
        self.assertEqual(tuple(red_tau.values.blocks), ((2, 3), (4, 3, 2, 1), (5, 2), (3,)))

        self.assertEqual(red_tau[1, 2], (2, 2))
        self.assertEqual(red_tau[0, 1], (4, 1))

    def test_pzero(self) -> None:
        # TODO
        pass

    def test_hash(self) -> None:
        G = LinearGroup((2, 3, 4))
        tau1 = Tau.from_flatten([6,2,1,4,1,2,5,3,1], G)
        tau2 = Tau.from_flatten([6,2,1,4,1,2,5,3,1], G)
        self.assertEqual(tau1, tau2)
        self.assertEqual(hash(tau1), hash(tau2))

        rtau1 = tau1.reduced
        rtau2 = tau2.reduced
        self.assertEqual(rtau1, rtau2)
        self.assertEqual(hash(rtau1), hash(rtau2))
