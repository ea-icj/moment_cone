import unittest
from typing import Any

from moment_cone.representation import *

class TestRepresentation(unittest.TestCase):
    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        from moment_cone.utils import getLogger
        import logging
        getLogger().setLevel(logging.ERROR)

    def test_instance_base(self) -> None:
        # Testing instance with same construction
        V1 = KroneckerRepresentation((5, 5, 5, 1), seed=0)
        V2 = KroneckerRepresentation((5, 5, 5, 1), seed=0)
        self.assertIs(V1, V2)
        
    def test_instance_Kro_dim(self) -> None:
        # Testing that it doesn't depend on the last dimension for Kronecker
        V1 = KroneckerRepresentation((5, 5, 5), seed=0)
        V2 = KroneckerRepresentation((5, 5, 5, 1), seed=0)
        self.assertIs(V1, V2)

    def test_instance_seed(self) -> None:
        # Testing that reusing seed allows to get back the instance
        V1 = KroneckerRepresentation((5, 5, 5, 1))
        V2 = KroneckerRepresentation((5, 5, 5, 1), seed=V1.seed)
        self.assertIs(V1, V2)

    def test_instance_optional(self) -> None:
        # Testing that it doesn't depend on the overridden parameters
        V1 = KroneckerRepresentation((5, 5, 5, 1), seed=0)
        V2 = KroneckerRepresentation((5, 5, 5, 1), seed=0, random_deep=1)
        self.assertIs(V1, V2)

    def test_instance_kro_particle(self) -> None:
        # Testing that Kronecker and Particle with same dimension doesn't
        # return same instance
        V1 = KroneckerRepresentation((5,), seed=0)
        V2 = FermionRepresentation((5,), particle_cnt=1, seed=0)
        self.assertIsNot(V1, V2)

    def test_instance_boson_fermion(self) -> None:
        # Testing that different Particle representation
        # leads to different instances
        V1 = BosonRepresentation((5,), particle_cnt=1, seed=0)
        V2 = FermionRepresentation((5,), particle_cnt=1, seed=0)
        self.assertIsNot(V1, V2)

    def test_seed(self) -> None:
        # Ensuring that seed differs for different instances
        V1 = KroneckerRepresentation((5, 5, 5, 1))
        V2 = KroneckerRepresentation((5, 5, 5, 1))
        self.assertIsNot(V1, V2)
        self.assertNotEqual(V1._generate_seed("unittest"), V2._generate_seed("unittest"))
        
    def test_seed_boson_fermion(self) -> None:
        # Ensuring that seed is different per particle representation
        V1 = BosonRepresentation((5,), particle_cnt=1, seed=0)
        V2 = FermionRepresentation((5,), particle_cnt=1, seed=0)
        self.assertNotEqual(V1._generate_seed("unittest"), V2._generate_seed("unittest"))

    def test_seed_particles(self) -> None:
        # Ensuring that seed is different when particle number change
        V1 = BosonRepresentation((5,), particle_cnt=1, seed=0)
        V2 = BosonRepresentation((5,), particle_cnt=2, seed=0)
        self.assertNotEqual(V1._generate_seed("unittest"), V2._generate_seed("unittest"))
