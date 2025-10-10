__all__ = (
    "Representation",
    "KroneckerRepresentation",
    "ParticleRepresentation",
    "BosonRepresentation",
    "FermionRepresentation",
)

from abc import ABC, abstractmethod
from functools import cached_property
import itertools
import numpy as np
from numpy.typing import NDArray
from typing import NamedTuple
from argparse import ArgumentParser, Namespace

from .typing import *
from .linear_group import LinearGroup
from .weight import Weight as WeightBase, WeightAsList, WeightAsListOfList
from .partition import Partition
from .rings import Matrix, Vector, Ring, PolynomialRingForWeights,PolynomialRing, Polynomial, Variable, I
from .root import Root


class TPi3DResult(NamedTuple):
    """ Result class of Representation.T_Pi_3D method """
    Q: NDArray[np.int64]
    QI: NDArray[np.int64]
    QV: NDArray[Any]
    QV_int: NDArray[Any]
    line_Q: NDArray[Any]
    line_QV: NDArray[Any]
    dict_Q: list[dict[Polynomial, Polynomial]]
    dict_QV: list[dict[Polynomial, Polynomial]] # FIXME: name and type is not correct anymore

    @overload
    def __call__(
            self,
            method: Method,
            kind: Optional[Literal["imaginary", "line"]] = None,
        ) -> NDArray[Any]:
        ...
    
    @overload
    def __call__(
            self,
            method: Method,
            kind: Literal["dict"]
        ) -> list[dict[Polynomial, Polynomial]]:
        ...

    def __call__(
            self,
            method: Method,
            kind: Optional[Literal["imaginary", "line", "dict"]] = None,
        ) -> NDArray[Any] | list[dict[Polynomial, Polynomial]]:
        match kind, method:
            case None, 'probabilistic': return self.Q
            case 'imaginary', 'probabilistic': return self.QI
            case None | 'imaginary', 'symbolic': return self.QV
            case None, 'symbolic_int': return self.QV_int
            case 'line', 'probabilistic': return self.line_Q
            case 'line', 'symbolic': return self.line_QV
            case 'dict', 'probabilistic': return self.dict_Q
            case 'dict', 'symbolic': return self.dict_QV
            case _: raise ValueError(f"Unknown method {method} and kind {kind}")

    
class Representation(ABC):
    """ Base class of a representation """
    #: Cache of instances of Representation
    __all_instances: ClassVar[
        dict[Any, Self]
    ] = {}

    Weight: ClassVar[type[WeightBase]] = WeightBase #: Weight class
    G: LinearGroup
    random_deep: int #: Deepness of the probabilistic methods
    seed: int #: Seed for the pseudo-random generators

    def __new__(cls: type[Self], 
                G: LinearGroup | Iterable[int],
                *,
                random_deep: int = 1,
                seed: Optional[int] = None,
                extra: Optional[Hashable] = None) -> Self:
        """ Custom instance construction in order to reuse previous instance """
        if not isinstance(G, LinearGroup):
            G = LinearGroup(G)

        # Generate specific seed
        from .utils import generate_seed
        seed = generate_seed(seed)

        # Get or create unique instance
        self_keys = (cls.__name__, tuple(G), seed, random_deep, extra)
        try:
            self = cls.__all_instances[self_keys]
        except KeyError:
            self = super().__new__(cls)
            self.G = G
            self.random_deep = random_deep
            self.seed = seed
            cls.__all_instances[self_keys] = self
        return self
    
    def _generate_seed(self, extra: Optional[Hashable] = None, nbytes: int = 8) -> int:
        """ Generate sub-seed specific to the representation
        
        >>> V1 = KroneckerRepresentation((3, 3, 3, 1), seed=888)
        >>> seed1 = V1._generate_seed("doctest")
        >>> seed2 = V1._generate_seed("other")
        >>> seed1 == seed2
        False

        >>> V2 = KroneckerRepresentation((3, 3, 3, 1), seed=777)
        >>> seed3 = V2._generate_seed("doctest")
        >>> seed1 == seed3
        False

        >>> V3 = KroneckerRepresentation((3, 3, 3, 1), seed=888)
        >>> seed4 = V3._generate_seed("doctest")
        >>> seed1 == seed4
        True
        """
        from .utils import generate_seed
        self_keys = (
            type(self).__name__,
            tuple(self.G),
            self.random_deep,
            extra
        )
        return generate_seed(self.seed, self_keys, nbytes)

    def _manual_seed(self, extra: Optional[Hashable] = None, nbytes: int = 8) -> None:
        """ Set the seed for used pseudo-random generators
        
        >>> import random
        >>> V = KroneckerRepresentation((3, 3, 3, 1))
        >>> V._manual_seed("doctest")
        >>> r1 = random.random()

        >>> V._manual_seed("other")
        >>> r2 = random.random()
        >>> r1 == r2
        False

        >>> V._manual_seed("doctest")
        >>> r3 = random.random()
        >>> r1 == r3
        True
        """
        from .utils import manual_seed
        seed = self._generate_seed(extra, nbytes)
        manual_seed(seed)

    def reduce(self, reduced_G: Iterable[int] | LinearGroup, **extra: Any) -> Self:
        """ Return a similar representation but with modified linear group

        Ensure that extra properties like seed, random_deep, etc. are
        set so that to ensure proper unique instance recovery.

        >>> V = KroneckerRepresentation((4, 4, 4, 1), random_deep=2)
        >>> Vred1 = V.reduce((4, 3, 3, 1))
        >>> Vred2 = V.reduce((4, 3, 3, 1))
        >>> Vred1 is Vred2
        True

        >>> V = FermionRepresentation((8,), particle_cnt=3, random_deep=2)
        >>> Vred1 = V.reduce((6,))
        >>> Vred2 = V.reduce((6,))
        >>> Vred1 is Vred2
        True

        >>> V = BosonRepresentation((8,), particle_cnt=3, random_deep=2)
        >>> Vred1 = V.reduce((6,), particle_cnt=2)
        >>> Vred2 = V.reduce((6,), particle_cnt=2)
        >>> Vred1.particle_cnt
        2
        >>> Vred1 is Vred2
        True
        """
        kwargs: dict[str, Any] = dict(
            random_deep=self.random_deep,
            seed=self.seed,
        )
        kwargs.update(extra)
        return type(self)(
            reduced_G,
            **kwargs,
        )
    
    def __getnewargs_ex__(self) -> tuple[tuple[Any, ...], dict[str, Any]]:
        """ Minimal state that need to be passed to __new__ in order to get the proper instance """
        return (), dict(
            G=self.G,
            random_deep=self.random_deep,
            seed=self.seed,
        )

    def __getstate__(self) -> dict[str, Any]:
        """ Minimal state that reproduce the instance (for serialization) """
        return dict(
            G=self.G,
            random_deep=self.random_deep,
            seed=self.seed,
        )
    
    def __setstate__(self, state: dict[str, Any]) -> None:
        """ Restoring instance from it's serialization state """
        self.G = state["G"]
        self.random_deep = state["random_deep"]
        self.seed = state["seed"]

    def weight(self, *args: Any, **kwargs: Any) -> WeightBase:
        """ Creates a weight for the given representation """
        return self.Weight(self.G, *args, **kwargs)
    
    @cached_property
    @abstractmethod
    def all_weights(self) -> Sequence[WeightBase]:
        """ All meaningful weights in a specific order """
        ...

    @abstractmethod
    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        """ Index of a given weight in the sequence returned by all_weights

        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """
        ...

    @cached_property
    @abstractmethod
    def dim_cone(self) -> int:
        """ Expected dimension of the cone. To be checked with Stabilizer of K """
        ...

    @cached_property
    @abstractmethod
    def dim(self) -> int:
        """ Dimension of V """
        ...

    @property
    @abstractmethod
    def weights_mod_outer(self) -> Iterable[WeightBase]:
        """
        Returns all decreasing weights modulo the symmetries of V that is of G.
        
        It assumes that G is arranged so that blocks are contiguous.
        """
        ...

    @cached_property
    @abstractmethod
    def T_Pi_3D(self) -> TPi3DResult:
        """
        The list of matrices rho_V(xi) for xi in the bases of K as a tridimensional np.array.
        The first entry are indexed by all_rootsK using the dictionary dict_rootK of the class LinearGroup.
        The other entries are indexed by self.all_Weights using self.index_of_weight(chi).

        The result of this property is cached and is computer from random numbers.
        If you need to get new random elements, call `clear_T_Pi_3D` method before.
        """
    ...

    def clear_T_Pi_3D(self) -> None:
        """ Clear cache of T_Pi_3D property """
        try:
            del self.__dict__["T_Pi_3D"]
        except KeyError:
            pass


    @cached_property
    @abstractmethod
    def actionK(self) -> NDArray[Any]:
        """
        The list of matrices rho_V(xi) for xi in the bases of K as a tridimensional np.array.
        The first entry are indexed by all_rootsK using index_in_all_of_K of the class LinearGroup.
        The other entries are indexed by self.all_Weights using self.index_of_weight(chi).
        """
        ... 
                
                

    @abstractmethod
    def Matrix_Graph(self, roots: Iterable[Root]) -> Matrix:
        """
        Return the matrix of the graph indexed by self.all_weights
        """
        ...
    
    def __repr__(self) -> str:
        return f"{type(self).__name__}(G={self.G})"


    @cached_property
    def Q(self) -> Ring:
        from .rings import QQ
        return QQ
    
    @cached_property
    def QI(self) -> Ring:
        from .rings import QQ, I
        return QQ[I]

    @cached_property
    def QZ(self) -> PolynomialRingForWeights:
        from .rings import PolynomialRingForWeights, QQ
        return PolynomialRingForWeights(QQ, "z")
    
    @cached_property
    def QV(self) -> PolynomialRingForWeights:
        from .rings import PolynomialRingForWeights, QQ
        from .weight import Weight
        return PolynomialRingForWeights(QQ, weights=self.all_weights)
    
    @cached_property
    def QV2(self) -> PolynomialRingForWeights:
        from .rings import PolynomialRingForWeights, QQ
        from .weight import Weight
        return PolynomialRingForWeights(QQ,
            weights=self.all_weights,
            seed=('va', 'vb'),
        )

    @cached_property
    def QIV(self) -> PolynomialRingForWeights:
        from .rings import PolynomialRingForWeights, QQ, I
        from .weight import Weight
        return PolynomialRingForWeights(
            QQ[I], 
            weights=self.all_weights,
            seed=('vr', 'vi')
        )
    
    @cached_property
    def QU_QV(self) -> PolynomialRingForWeights:
        return (self.G).QU((self.QV.sage_ring).fraction_field())

    def random_element(self) -> NDArray[np.int64]:
        """ Random vector avoiding 0 entries """
        return (-1)**np.random.randint(0,2,size=self.dim)*np.random.randint(1, 10000, size=self.dim)
    
    @cached_property
    def fixed_random_element_Q(self) -> NDArray[np.int64]:
        from .utils import manual_seed
        self._manual_seed("fixed_random_element_Q")
        return self.random_element()
    
    @cached_property
    def fixed_random_element_QI(self) -> Polynomial:
        from .rings import I
        from .utils import manual_seed
        self._manual_seed("fixed_random_element_QI")
        return self.random_element() + I * self.random_element()

    @cached_property
    def fixed_random_line_in(self) -> Polynomial:
        from .utils import manual_seed
        self._manual_seed("fixed_random_line_in")
        return self.random_element() * self.QZ('z') + self.random_element()
    
    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments that defines the representation """
        from .utils import to_literal
        group = parent_parser.add_argument_group(
            "Representation"
        )
        group.add_argument(
            "representation",
            type=lambda s: to_literal(Literal["Kronecker", "Fermion", "Boson"], s),
            choices=("Kronecker", "Fermion", "Boson"),
            help="Representation type",
        )
        group.add_argument(
            "N",
            type=int,
            nargs='+',
            help="Dimensions of the linear groups",
        )
        group.add_argument(
            "--particle_cnt",
            type=int,
            default=None,
            help="Number of particles in Fermion and Boson representation",
        )
        group.add_argument(
            "--seed",
            type=int,
            default=None,
            help="Seed for the pseudo-random generators",
        )    
        group.add_argument(
            "--random_deep",
            type=int,
            default=1,
            help="Deepness of the check for probabilistic methods",
        )

    @classmethod
    def from_config(cls: type[Self], config: Namespace, **kwargs: Any) -> "Representation":
        """ Build a representation from the command-line arguments """
        G = LinearGroup(config.N)
        V: Representation

        common_args = dict(
            G=G,
            random_deep=config.random_deep,
            seed=config.seed,
        )
        specific_args: dict[str, Any] = {}
        repr_class: type[Representation]

        match config.representation.lower():
            case "kronecker":
                from .representation import KroneckerRepresentation
                repr_class = KroneckerRepresentation
            case "fermion":
                assert config.particle_cnt is not None, "particle_cnt is mandatory for Fermion representation"
                from .representation import FermionRepresentation
                repr_class = FermionRepresentation
                specific_args["particle_cnt"] = config.particle_cnt
            case "boson":
                assert config.particle_cnt is not None, "particle_cnt is mandatory for Boson representation"
                from .representation import BosonRepresentation
                repr_class = BosonRepresentation
                specific_args["particle_cnt"] = config.particle_cnt
            case _:
                raise ValueError(f"Invalid representation name {config.representation}")
        
        V = repr_class(**common_args, **specific_args)
        return V
    

class KroneckerRepresentation(Representation):
    Weight = WeightAsList
    
    def __new__(cls: type[Self], G: LinearGroup | Iterable[int], **kwargs: Any) -> Self:
        """ Kronecker representation

        The last dimension of the linear group must be one for consistency reason.

        >>> KroneckerRepresentation((3, 3, 3, 1))
        KroneckerRepresentation(G=GL(3)xGL(3)xGL(3)xGL(1))
        >>> KroneckerRepresentation((3, 3, 3))
        KroneckerRepresentation(G=GL(3)xGL(3)xGL(3)xGL(1))
        """
        G = LinearGroup(G)
        if G[-1] != 1:
            from .utils import getLogger
            logger = getLogger("KroneckerRepresentation")
            logger.warning("Dimension 1 appended to the linear group of the Kronecker representation")
            G = LinearGroup(tuple(G) + (1,))
        return super().__new__(cls, G, **kwargs)

    @cached_property
    def dim_cone(self) -> int:
        return self.G.rank - len(self.G) + 1
    
    @cached_property
    def dim(self) -> int:
        from .utils import prod
        return prod(self.G)
    
    @cached_property
    def all_weights(self) -> list[WeightAsList]:
        """
        Returns all possible weights in the lexicographical order
        """
        L: list[WeightAsList] = []
        for idx, w in enumerate(itertools.product(*(range(di) for di in self.G))):
            L.append(WeightAsList(self.G, as_list=w, index=idx))
        return L

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        """
        Returns index of the weight in the lexicographical order for the current representation (see `all_weights` method).
        
        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """
        if not use_internal_index or chi.index is None:
            if not isinstance(chi, WeightAsList):
                raise ValueError("Invalid weight representation")        
            from operator import mul
            stride = itertools.accumulate(reversed(self.G[1:]), mul, initial=1)
            chi.index = sum(v * s for v, s in zip(reversed(chi.as_list), stride))
        return chi.index

    @property
    def weights_mod_outer(self) -> Iterable[WeightAsList]:
        """
        Returns all decreasing weights modulo the symmetries of V that is of G.
        
        It assumes that G is arranged so that blocks are contiguous.
        """
        # Using Partition to generate all decreasing weights within a block.
        # The whole weights will be defined as the Cartesian product of the weights for each block.
        from .utils import group_by_block
        from .rings import vector
        
        def pad(p: Partition, l: int) -> tuple[int, ...]:
            return p.pad(l)
        
        block_weights = tuple(
            tuple( # Converting to tuple (instead of a map) seems necessary to keep the right hi (FIXME)
                p.pad(hi) # Adding trailing zeros if necessary
                for p in Partition.all_of_length(hi, di - 1)
            )
            for di, hi in group_by_block(self.G) # Compress returns (value, multiplicity) for each block of G
        )

        for w in itertools.product(*block_weights):
            Lw=[]
            for x in w :
                Lw+=list(x)
            yield WeightAsList(
                self.G,
                as_list=list(vector(sum(w, start=())))
            ) # Summing tuples is concatenating them 
    
    @cached_property
    def T_Pi_3D(self) -> TPi3DResult:
        """
        The list of matrices rho_V(xi) for xi in the bases of K as a tridimensional np.array.
        The first entry are indexed by all_rootsK using the dictionary dict_rootK of the class LinearGroup.
        The other entries are indexed by self.all_Weights using self.index_of_weight(chi).
        """
        from .utils import manual_seed
        self._manual_seed("T_Pi_3D")

        # Computation made once
        result_Q = np.zeros((self.random_deep,self.dim, self.dim, self.G.dimU), dtype=np.int64)
        result_QI = np.zeros((2*self.random_deep,self.dim, self.dim, self.G.dimU), dtype=np.int64) #first index is used for real and imaginary part.
        result_QV = np.zeros((self.dim, self.dim, self.G.dimU), dtype=object)
        result_QV_int = np.zeros((self.dim, self.dim, self.G.dimU), dtype=np.int8)
        result_line_Q = np.zeros((2*self.random_deep,self.dim, self.dim, self.G.dimU), dtype=np.int64)
        result_line_QV = np.zeros((self.dim, self.dim, self.G.dimU), dtype=object)
        
        K=self.QV2.fraction_field()
        ring_R0 = PolynomialRing(K,"z")
        dict_Q: list[dict[Polynomial, Polynomial]] = [{} for _ in range(self.random_deep)]
        #dict_QV: list[dict[Polynomial, Polynomial]] = [{} for _ in range(self.random_deep)]
        dict_QV: list[dict[Polynomial, Polynomial]] = [{}]
        # produce a collection of 5* random_deep random vectors 
        random_vectors =(-1)**np.random.randint(0,2,size=(5*self.random_deep,self.dim))*np.random.randint(1, 1000, size=(5*self.random_deep,self.dim))
        # Index 0 used for Q, 1 and 2, for QI (real and imaginary parts), 3,4 for line_Q (a and b for az+b)

        for chi in self.all_weights:
            id_chi=self.index_of_weight(chi)
            vchi_a, vchi_b = self.QV2.variable(chi)
            #dict_QV[0][self.QV.variable(chi)]= vchi_a*ring_R0('z') + vchi_b # type: ignore
            #homs_QV[0][self.QV.variable(chi)]= vchi_a*ring_R0('z') + vchi_b # type: ignore
            #_QV.append(vchi_a*ring_R0('z') + vchi_b)#coucou
            dict_QV[0][self.QV.variable(chi)]=vchi_a*ring_R0('z') + vchi_b # type: ignore
            for k,b in enumerate(chi.as_list):
                for i in range(b):
                    chi_i = WeightAsList(
                        self.G,
                        as_list=chi.as_list[:k] + (i,) + chi.as_list[k+1:]
                        )
                    id_i = self.index_of_weight(chi_i)
                    for p in range(self.random_deep):
                        result_Q[p,id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = random_vectors[5*p,id_chi]
                        result_QI[2*p,id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = random_vectors[5*p+1,id_chi]
                        result_QI[2*p+1,id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = random_vectors[5*p+2,id_chi]
                        result_QV[id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = self.QV.variable(chi)
                        result_QV_int[id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = 1
                        result_line_Q[2*p,id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = random_vectors[5*p+3,id_chi]
                        result_line_Q[2*p+1,id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = random_vectors[5*p+4,id_chi]
                        dict_Q[p][self.QV.variable(chi)]= random_vectors[5*p+3,id_chi]*self.QZ('z')+random_vectors[5*p+4,id_chi]
                        result_line_QV[id_chi,id_i,Root(k,i,b).index_in_all_of_U(self.G)] = vchi_a*ring_R0('z') + vchi_b

        homs_Q=[]
        for p in range(self.random_deep):
            subs_Q = [dict_Q[p].get(self.QV.variable(chi),1) for chi in self.all_weights]
            homs_Q.append(self.QV.hom(subs_Q,self.QZ.sage_ring))
        #homs_QV=[self.QV.hom(subs_Q,self.QV2.sage_ring)] # List of length one since in deterministic case non need of repetition    
        return TPi3DResult(
            result_Q, result_QI, result_QV, result_QV_int,
            result_line_Q, result_line_QV,
            homs_Q, dict_QV
        )
    
    
    @cached_property
    def actionK(self) -> NDArray[Any]:
        """
        The list of matrices rho_V(xi) for xi in the bases of K as a tridimensional np.array.
        The first entry are indexed by all_rootsK using index_in_all_of_K of the class LinearGroup.
        The other entries are indexed by self.all_Weights using self.index_of_weight(chi).
        """
       
        shiftI = self.dim # basis over the real e_0,...,e_{D-1},Ie_0,Ie_1,...
        result=np.zeros((self.G.dim,2*self.dim,2*self.dim), dtype=np.int64)
        for chi in self.all_weights:
            id_chi=self.index_of_weight(chi)
            for k,b in enumerate(chi.as_list):
                # entries for action of I E^k_bb
                result[Root(k,b,b).index_in_all_of_K(self.G),shiftI+id_chi,id_chi]=1
                result[Root(k,b,b).index_in_all_of_K(self.G),id_chi,shiftI+id_chi]=-1
                for j in range(b+1,self.G[k]):
                    chi_j = WeightAsList(
                        self.G,
                        as_list=chi.as_list[:k] + (j,) + chi.as_list[k+1:]
                        )
                    id_j = self.index_of_weight(chi_j)
                    
                    # entries for action of E^k_ij - E^k_ji
                    result[Root(k,b,j).index_in_all_of_K(self.G),id_j,id_chi]=-1
                    result[Root(k,b,j).index_in_all_of_K(self.G),shiftI+id_j,shiftI+id_chi]=-1
                    result[Root(k,j,b).index_in_all_of_K(self.G),shiftI+id_j,id_chi]=1
                    result[Root(k,j,b).index_in_all_of_K(self.G),id_j,shiftI+id_chi]=-1

                for i in range(b):
                    chi_i = WeightAsList(
                        self.G,
                        as_list=chi.as_list[:k] + (i,) + chi.as_list[k+1:]
                        )
                    id_i = self.index_of_weight(chi_i)
                    result[Root(k,i,b).index_in_all_of_K(self.G),id_i,id_chi]=1
                    result[Root(k,i,b).index_in_all_of_K(self.G),shiftI+id_i,shiftI+id_chi]=1
                    result[Root(k,b,i).index_in_all_of_K(self.G),shiftI+id_i,id_chi]=1
                    result[Root(k,b,i).index_in_all_of_K(self.G),id_i,shiftI+id_chi]=-1
        return(result)            
    

    def Matrix_Graph(self, roots : Iterable[Root]) -> Matrix:
        """
        Return the matrix of the graph indexed by self.all_weights
        """
        from .rings import matrix, ZZ

        M = matrix(ZZ,self.dim,self.dim)
        for alpha in roots:
            Gred = LinearGroup(self.G[:alpha.k ] +self.G[alpha.k+1:])
            Vred = self.reduce(Gred)
            for w in Vred.all_weights:
                wj = WeightAsList(
                    self.G,
                    as_list=w.as_list[:alpha.k] + (alpha.j,) + w.as_list[alpha.k:]
                )
                idj = self.index_of_weight(wj)                   
                wi = WeightAsList(
                    self.G,
                    as_list=w.as_list[:alpha.k] + (alpha.i,) + w.as_list[alpha.k:]
                )
                idi = self.index_of_weight(wi)
                M[idi, idj]=1
        return M


class ParticleRepresentation(Representation):
    """ Representation specific to physical particles """
    particle_cnt: int

    def __new__(
            cls: type[Self],
            G: LinearGroup | Iterable[int],
            *,
            particle_cnt: int,
            **kwargs: Any) -> Self:
        G = LinearGroup(G)
        if len(G) != 1:
            raise NotImplementedError("Product of GL not supported for particle representation")
        self = super().__new__(cls, G, **kwargs, extra=particle_cnt)
        self.particle_cnt = particle_cnt
        return self

    def _generate_seed(self, extra: Optional[Hashable] = None, nbytes: int = 8) -> int:
        return super()._generate_seed(
            extra=(self.particle_cnt, extra),
            nbytes=nbytes
        )

    def reduce(self, reduced_G: Iterable[int] | LinearGroup, **extra: Any) -> Self:
        kwargs: dict[str, Any] = dict(particle_cnt=self.particle_cnt)
        kwargs.update(extra)
        return super().reduce(reduced_G, **kwargs)
        
    def __getnewargs_ex__(self) -> tuple[tuple[Any, ...], dict[str, Any]]:
        """ Minimal state that need to be passed to __new__ in order to get the proper instance """
        args, kwargs = super().__getnewargs_ex__()
        kwargs["particle_cnt"] = self.particle_cnt
        return args, kwargs
    
    def __getstate__(self) -> dict[str, Any]:
        """ Minimal state that reproduce the instance (for serialization) """
        state = super().__getstate__()
        state["particle_cnt"] = self.particle_cnt
        return state
    
    def __setstate__(self, state: dict[str, Any]) -> None:
        """ Restoring instance from it's serialization state """
        super().__setstate__(state)
        self.particle_cnt = state["particle_cnt"]

    @cached_property
    def dim_cone(self) -> int:
        return self.G.rank

    @cached_property
    @abstractmethod
    def all_weights(self) -> Sequence[WeightAsListOfList]:
        ...

    @property
    def weights_mod_outer(self) -> Iterable[WeightAsListOfList]:
        return self.all_weights
        
    def weights_of_S(self, p: Partition) -> Iterable[WeightAsListOfList] : # Could be improved
        """
        Create the list of weights of S subset T acting on V. With multiplicities.
        S is given by p.
        Only used for V.type='fermion' or 'boson'
        """
        from math import comb
        from .utils import unique_combinations, prod

        GS = LinearGroup([len(p)])
        List_indices = []
        for i,x in enumerate(p) :
            List_indices += x*[i]
        
        chiTot = WeightAsListOfList(GS, as_list_of_list=[List_indices])
        for i,l in enumerate(unique_combinations(List_indices, self.particle_cnt)):
            chi = WeightAsListOfList(GS, as_list_of_list=[l], index=i)
            chi.mult = prod([
                comb(chiTot.as_vector[i], chi.as_vector[i])
                for i in range(len(p))
            ])
            yield chi

    @cached_property
    def actionK(self) -> NDArray[Any]:
        """
        The list of matrices rho_V(xi) for xi in the bases of K as a tridimensional np.array.
        The first entry are indexed by all_rootsK using index_in_all_of_K of the class LinearGroup.
        The other entries are indexed by self.all_Weights using self.index_of_weight(chi).
        """
        
        shiftI = self.dim # basis over the real e_0,...,e_{D-1},Ie_0,Ie_1,...
        result=np.zeros((self.G.dim,2*self.dim,2*self.dim), dtype=np.int64)
        
        for chi in self.all_weights:
            id_chi=self.index_of_weight(chi)
            for k,b in enumerate(chi.as_list_of_list[0]):
                index_b = chi.as_list_of_list[0].index(b)
                if k == index_b : # otherwise already done
                    mult = chi.as_list_of_list[0].count(b)
                    # entries for action of I E^k_bb
                    result[Root(0,b,b).index_in_all_of_K(self.G),shiftI+id_chi,id_chi]=1
                    result[Root(0,b,b).index_in_all_of_K(self.G),id_chi,shiftI+id_chi]=-1
                
                    # split chi 
                    L1 = chi.as_list_of_list[0][:index_b]
                    L2 = chi.as_list_of_list[0][index_b+1:]
                    for j in range(b+1, self.G[0]):
                        if isinstance(self, BosonRepresentation) or j not in L2:
                            L3 = tuple(sorted((j,) + L2))
                            if isinstance(self, BosonRepresentation) :
                                dec=0
                            else :
                                dec=L3.index(j)
                            Lj=L1+L3
                            chi_j = WeightAsListOfList(self.G, as_list_of_list=[Lj])
                            id_j = self.index_of_weight(chi_j)
                            # entries for action of E^k_ij - E^k_ji
                            result[Root(0,b,j).index_in_all_of_K(self.G),id_j,id_chi]=-mult* (-1)**dec
                            result[Root(0,b,j).index_in_all_of_K(self.G),shiftI+id_j,shiftI+id_chi]=-mult* (-1)**dec
                            result[Root(0,j,b).index_in_all_of_K(self.G),shiftI+id_j,id_chi]=mult* (-1)**dec
                            result[Root(0,j,b).index_in_all_of_K(self.G),id_j,shiftI+id_chi]=-mult* (-1)**dec

                    for i in range(b):
                        if isinstance(self, BosonRepresentation) or i not in L1:
                            L3 = tuple(sorted(L1 + (i,)))
                            if isinstance(self, BosonRepresentation) :
                                dec=0
                            else :
                                dec=len(L3)-L3.index(i)-1
                            Li=L3+L2
                            chi_i = WeightAsListOfList(self.G, as_list_of_list=[Li])
                            id_i = self.index_of_weight(chi_i) 
                            result[Root(0,i,b).index_in_all_of_K(self.G), id_i, id_chi] = mult* (-1)**dec
                            result[Root(0,i,b).index_in_all_of_K(self.G), shiftI + id_i, shiftI + id_chi] = mult* (-1)**dec
                            result[Root(0,b,i).index_in_all_of_K(self.G), shiftI + id_i, id_chi] = mult* (-1)**dec
                            result[Root(0,b,i).index_in_all_of_K(self.G), id_i, shiftI + id_chi] = -mult* (-1)**dec
        return(result)                    

    @cached_property
    def T_Pi_3D(self) -> TPi3DResult:
        """
        The list of matrices rho_V(xi) for xi in the bases of K as a tridimensional np.array.
        The first entry are indexed by all_rootsK using the dictionary dict_rootK of the class LinearGroup.
        The other entries are indexed by self.all_Weights using self.index_of_weight(chi).
        """
        from .utils import manual_seed
        self._manual_seed("T_Pi_3D")

        # Computation made once
        result_Q = np.zeros((self.random_deep,self.dim, self.dim, self.G.dimU), dtype=np.int64)
        result_QI = np.zeros((2*self.random_deep,self.dim, self.dim, self.G.dimU), np.int64)
        result_QV = np.zeros((self.dim, self.dim, self.G.dimU), dtype=object)
        result_QV_int = np.zeros((self.dim, self.dim, self.G.dimU), dtype=np.int8)
        result_line_Q = np.zeros((2*self.random_deep,self.dim, self.dim, self.G.dimU), dtype=np.int16)
        result_line_QV = np.zeros((self.dim, self.dim, self.G.dimU), dtype=object)        
        K=self.QV2.fraction_field()
        ring_R0 = PolynomialRing(K,"z")
        dict_Q: list[dict[Polynomial, Polynomial]] = [{} for _ in range(self.random_deep)]
        dict_QV=[{}]  # type: ignore 
        
        
        # produce a collection of 5* random_deep random vectors 
        random_vectors =(-1)**np.random.randint(0,2,size=(5*self.random_deep,self.dim))*np.random.randint(1, 1000, size=(5*self.random_deep,self.dim))
        # Index 0 used for Q, 1 and 2, for QI (real and imaginary parts), 3,4 for line_Q (a and b for az+b)

        for chi in self.all_weights:
            id_chi=self.index_of_weight(chi)
            vchi_a, vchi_b = self.QV2.variable(chi)
            dict_QV[0][self.QV.variable(chi)]=vchi_a*ring_R0('z') + vchi_b
            for k,b in enumerate(chi.as_list_of_list[0]):
                index_b = chi.as_list_of_list[0].index(b) #Used to treat repritions in the bosonic case
                if k == index_b:
                    mult = chi.as_list_of_list[0].count(b) #Constant obtained by derivative
                    # ssplit chi 
                    L1 = chi.as_list_of_list[0][:index_b]
                    L2 = chi.as_list_of_list[0][index_b+1:]
                    for i in range(b):
                        if isinstance(self, BosonRepresentation) or i not in L1: # otherwise action is zero
                            L3 = tuple(sorted(L1 + (i,)))
                            # dec is the number of moves used to resort i. Induces a sign in fermionic case.
                            if isinstance(self, BosonRepresentation) :
                                dec=0
                            else :
                                dec=len(L3)-L3.index(i)-1
                            Li=L3+L2
                            chi_i = WeightAsListOfList(self.G, as_list_of_list=[Li])
                            id_i = self.index_of_weight(chi_i)
                            for p in range(self.random_deep):
                                result_Q[p,id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult* (-1)**dec*random_vectors[5*p,id_chi]
                                result_QI[2*p,id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult*(-1)**dec*random_vectors[5*p+1,id_chi]
                                result_QI[2*p+1,id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult*(-1)**dec*random_vectors[5*p+2,id_chi]
                                result_QV[id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult* (-1)**dec*self.QV.variable(chi)
                                result_QV_int[id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult* (-1)**dec
                                result_line_Q[2*p,id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult* (-1)**dec*random_vectors[5*p+3,id_chi]  #(va[id_chi]*self.QZ('z')+vb[id_chi])
                                result_line_Q[2*p+1,id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult* (-1)**dec*random_vectors[5*p+4,id_chi]
                                dict_Q[p][self.QV.variable(chi)] = random_vectors[5*p+3,id_chi]*self.QZ('z')+random_vectors[5*p+4,id_chi]
                                result_line_QV[id_chi,id_i,Root(0,i,b).index_in_all_of_U(self.G)] = mult* (-1)**dec*(vchi_a*ring_R0('z') + vchi_b)
        homs_Q=[]
        for p in range(self.random_deep):
            subs_Q = [dict_Q[p].get(self.QV.variable(chi),1) for chi in self.all_weights]
            homs_Q.append(self.QV.hom(subs_Q,self.QZ.sage_ring))
        #homs_QV=[self.QV.hom(subs_Q,self.QV2.sage_ring)]    
        return TPi3DResult(
            result_Q, result_QI, result_QV,result_QV_int,
            result_line_Q, result_line_QV,
            homs_Q, dict_QV
        )
      

    def Matrix_Graph(self, roots : Iterable[Root]) -> Matrix:
        """
        Return the matrix of the graph indexed by self.all_weights
        """
        from .rings import matrix, ZZ

        M = matrix(ZZ,self.dim,self.dim)

        if isinstance(self, FermionRepresentation):
            shiftrank = 1
        else :
            shiftrank = 0

        Vred = self.reduce(
            LinearGroup([self.G[0] - shiftrank]),
            particle_cnt=self.particle_cnt - 1,
        )

        for alpha in roots:
            for w in Vred.all_weights: 
                L1=[s for s in w.as_list_of_list[0] if s<alpha.j]                
                L2=[s+shiftrank for s in w.as_list_of_list[0] if s>=alpha.j]
                lj=L1+[alpha.j]+L2  # we insert j
                if isinstance(self, BosonRepresentation) or alpha.i not in lj: # Otherwise E_ij v =0
                    wj  = WeightAsListOfList(self.G,as_list_of_list=[lj])
                    idj = self.index_of_weight(wj)
                    li=L1+[alpha.i]+L2  # we insert i
                    li.sort()
                    wi = WeightAsListOfList(self.G,as_list_of_list=[li])
                    idi = self.index_of_weight(wi)
                    if isinstance(self, FermionRepresentation):
                        M[idi, idj] = 1
                    else :
                        M[idi, idj] = 1
        return M

    def __repr__(self) -> str:
        return f"{type(self).__name__}(G={self.G}, particle_cnt={self.particle_cnt})"


class FermionRepresentation(ParticleRepresentation):
    """
    Representation specific to Fermion particles
    
    TODO: doctest
    """
    Weight = WeightAsListOfList

    @cached_property
    def dim(self) -> int:
        from math import comb
        return comb(self.G.rank, self.particle_cnt)

    @cached_property
    def all_weights(self) -> list[WeightAsListOfList]:
        L: list[WeightAsListOfList] = []
        weights = itertools.combinations(range(self.G[0]), self.particle_cnt)
        for i, w in enumerate(weights):
            L.append(WeightAsListOfList(self.G, as_list_of_list=(w,), index=i))
        return L

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        if not isinstance(chi, WeightAsListOfList):
            raise ValueError("Invalid weight representation")
        
        if not use_internal_index or chi.index is None:
            from math import comb
            S = chi.as_list_of_list[0]
            k = self.particle_cnt
            id=sum([comb(self.G.rank-j-1,k-1) for j in range(S[0])])# Subset with smaller first element
            for i,p in enumerate(itertools.pairwise(S)):
                id+=sum([comb(self.G.rank-j-1,k-i-2) for j in range(p[0]+1,p[1])]) # With j in position i+1 and equal for smaller indices       
            chi.index = id

        return chi.index


class BosonRepresentation(ParticleRepresentation):
    """ Representation specific to Boson particles """
    Weight = WeightAsListOfList

    @cached_property
    def dim(self) -> int:
        from math import comb
        return comb(
            self.G.rank + self.particle_cnt - 1,
            self.G.rank - 1
        )
    
    @cached_property
    def all_weights(self) -> list[WeightAsListOfList]:
        from .rings import vector, ZZ
        L: list[WeightAsListOfList] = []
        comb = itertools.combinations(
            range(self.G[0] + self.particle_cnt-1),
            self.G[0] - 1
        )
        for i, w in enumerate(comb):
            v = vector(ZZ, self.G.rank)
            v[0] = w[0]
            AsL = w[0]*[0]
            for j, p in enumerate(itertools.pairwise(w)):
                v[j + 1] = p[1] - p[0] - 1
                AsL += v[j + 1] * [j + 1]
            v[-1] = self.G.rank + self.particle_cnt - w[-1] - 2
            AsL += v[-1] * [self.G.rank - 1]
            L.append(WeightAsListOfList(
                self.G,
                as_vector=v,
                as_list_of_list=[AsL],
                index=i
            ))
        return L

    def index_of_weight(self, chi: WeightBase, use_internal_index: bool = True) -> int:
        if not isinstance(chi, WeightAsListOfList):
            raise ValueError("Invalid weight representation")
        
        if not use_internal_index or chi.index is None:
            for i, chi2 in enumerate(self.all_weights):
                if chi2.as_vector == chi.as_vector :
                    chi.index = i
                    break
            assert chi.index is not None
            
        return chi.index
