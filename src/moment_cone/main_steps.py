
__all__ = (
    "Dataset",
    "ListDataset",
    "Step",
    "GeneratorStep",
    "FilterStep",
    "TransformerStep",
    "GeneralStabilizerDimensionCheck",
    "TauCandidatesStep",
    "SubModuleConditionStep",
    "StabilizerConditionStep",
    "InequalityCandidatesStep",
    "PiDominancyStep",
    "LinearTriangularStep",
    "BKRConditionStep",
    "BirationalityStep",
    "GrobnerStep",
    "ExportStep",
    "MomentConeStep",
)

from argparse import ArgumentParser, Namespace
import typing
from tqdm.auto import tqdm

from .typing import *
from .representation import Representation
from .linear_group import LinearGroup
from .tau import Tau
from .inequality import Inequality
from .root import Root
from .kronecker import KroneckerCoefficient, KroneckerCoefficientMLCache
from .bkr import PlethysmCache
from .utils import to_literal
from .export import ExportFormat
from .tau_storage import UniqueTauStr

class Dataset(Generic[T], ABC):
    """ Catalog of pending and validated objects of type T
    
    T will be typically Tau or Inequality

    Pending objects (Tau or Inequality) are those that are not excluded from
    the process but still need to be validated by a filter.

    Validated objects are those that we are sure to be valid
    and thus should always be keeped by the filters.
    """
    @abstractmethod
    def __init__(
            self,
            pending_or_all: Iterable[T] | Iterable[tuple[T, bool]], #: pending elements or joined pending/validated elements as tuple with status
            validated: Optional[Iterable[T]] = None #: validated elements or None if the first elements are tuple of elements and status
        ):
        """
        Two construction ways:
        1. giving pending and validated elements separately using both arguments
        2. given pending & validated elements as tuple (element, status) by giving only the first parameter
        
        Status: True for validated, False for pending.
        """
        ...

    @classmethod
    def from_separate(cls, pending: Iterable[T] = [], validated: Iterable[T] = []) -> Self:
        return cls(pending, validated)
    
    @classmethod
    def from_all(cls, all_elements: Iterable[tuple[T, bool]] = []) -> Self:
        return cls(all_elements)

    @abstractmethod
    def pending(self) -> Iterable[T]:
        ...
    
    @abstractmethod
    def validated(self) -> Iterable[T]:
        ...

    def all(self) -> Iterable[tuple[T, bool]]:
        """ All elements as tuple (element, status)
        
        Status: True for validated, False for pending.
        """
        from itertools import chain
        return chain(
            map(lambda t: (t, False), self.pending()),
            map(lambda t: (t, True), self.validated())
        )

    def __iter__(self) -> Iterator[T]:
        from itertools import chain
        return chain(self.pending(), self.validated())
    
    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


class ListDataset(Dataset[T]):
    """ Catalog of pending/validated objects implemented using list """
    __pending: list[T]
    __validated: list[T]

    def __init__(
            self,
            pending_or_all: Iterable[T] | Iterable[tuple[T, bool]], #: pending elements or joined pending/validated elements as tuple with status
            validated: Optional[Iterable[T]] = None #: validated elements or None if the first elements are tuple of elements and status
        ):
        if validated is None:
           self.__pending = []
           self.__validated = []
           for element, status in cast(Iterable[tuple[T, bool]], pending_or_all):
               if status:
                   self.__validated.append(element)
               else:
                   self.__pending.append(element)
        else:
            pending = cast(Iterable[T], pending_or_all)
            self.__pending = pending if isinstance(pending, list) else list(pending)
            self.__validated = validated if isinstance(validated, list) else list(validated)

    def pending(self) -> list[T]:
        return self.__pending
    
    def validated(self) -> list[T]:
        return self.__validated

    def __repr__(self) -> str:
        return f"ListDataset(#pending={len(self.pending())}, #validated={len(self.validated())})"


class LazyDataset(Dataset[T]):
    """ Catalog of pending/validated objects implemented in a lazy way

    That means that elements are stored only if necessary (buffer) and that this
    dataset cannot be iterated twice.
    """
    __all: Iterator[tuple[T, bool]]
    __pending_cache: list[T]
    __validated_cache: list[T]
    has_ended: bool
    pending_cnt: int
    validated_cnt: int

    def __init__(
            self,
            pending_or_all: Iterable[T] | Iterable[tuple[T, bool]], #: pending elements or joined pending/validated elements as tuple with status
            validated: Optional[Iterable[T]] = None #: validated elements or None if the first elements are tuple of elements and status
        ):
        """
        Two construction ways:
        1. giving pending and validated elements separately using both arguments
        2. given pending & validated elements as tuple (element, status) by giving only the first parameter
        
        Status: True for validated, False for pending.
        """
        if validated is None:
            self.__all = iter(map(
                self.__count_status,
                cast(Iterable[tuple[T, bool]], pending_or_all
            )))
        else:
            from itertools import chain
            self.__all = iter(map(self.__count_status,
                chain(
                    map(lambda t: (t, False), cast(Iterable[T], pending_or_all)),
                    map(lambda t: (t, True), validated)
                )
            ))

        self.__pending_cache = []
        self.__validated_cache = []
        self.pending_cnt = 0
        self.validated_cnt = 0
        self.has_ended = False

    def __count_status(self, el_and_status: tuple[T, bool]) -> tuple[T, bool]:
        element, status = el_and_status
        if status:
            self.validated_cnt += 1
        else:
            self.pending_cnt += 1
        return el_and_status

    def all(self) -> Iterable[tuple[T, bool]]:
        return self.__all
    
    def pending(self) -> Iterable[T]:
        while True:
            if len(self.__pending_cache):
                yield self.__pending_cache.pop()
            else:
                try:
                    t, validated = next(self.__all)
                except StopIteration:
                    self.has_ended = True
                    return
                if validated:
                    self.__validated_cache.append(t)
                else:
                    yield t
    
    def validated(self) -> Iterable[T]:
        while True:
            if len(self.__validated_cache):
                yield self.__validated_cache.pop()
            else:
                try:
                    t, validated = next(self.__all)
                except StopIteration:
                    self.has_ended = True
                    return
                if validated:
                    yield t
                else:
                    self.__pending_cache.append(t)

    def __repr__(self) -> str:
        continuing = "" if self.has_ended else "..."
        return f"LazyDataset(#pending={self.pending_cnt}{continuing}, #validated={self.validated_cnt}{continuing})"


class Step:
    """ Represents one computational step (generation, filtering, etc)
    
    A step is associated to a representation.

    It can add arguments to command-line interface and be construct from
    the parsed command-line arguments.

    Processing is done through __call__ (overload of () operator).
    """
    V: Representation
    TDataset: type[Dataset[Any]]

    def __init__(
            self,
            V: Representation,
            dataset_type: type[Dataset[Any]] = ListDataset,
            quiet: bool = False,
            **kwargs: Any
        ):
        self.V = V
        self.quiet = quiet
        self.TDataset = dataset_type

    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        pass

    @classmethod
    def from_config(
            cls: type[Self],
            V: Representation,
            config: Namespace,
            **kwargs: Any,
        ) -> Self:
        """ Build a step from the representation and the command-line arguments """
        try:
            return cls(V, config=config, quiet=config.quiet, **kwargs)
        except AttributeError:
            return cls(V, config=config, **kwargs)

    @property
    def name(self) -> str:
        return type(self).__name__
    
    @property
    def G(self) -> LinearGroup:
        return self.V.G
    
    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        """ Effective computation of the step """
        from tqdm.contrib.logging import logging_redirect_tqdm
        with logging_redirect_tqdm():
            self.apply(*args, **kwargs)

    @abstractmethod
    def apply(self, *args: Any, **kwargs: Any) -> Any:
        ...

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"
    
    def _tqdm(self,
              iterable: Optional[Iterable[T]] = None,
              desc: Optional[str] = None,
              leave: bool = False,
              disable: Optional[bool] = None,
              unit: str = 'it',
              **kwargs: Any,
              ) -> tqdm: # type: ignore
        """ Helper function to generate a progress bar with appropriate configuration """
        if desc is None:
            desc = type(self).__name__
        if disable is None:
            disable = self.quiet

        if (not disable or not self.quiet) and iterable is not None:
            from .utils import getLogger, IterableHook
            import logging
            logger = getLogger(type(self).__name__)
            if logger.isEnabledFor(logging.DEBUG):
                iterable = IterableHook(
                    iterable,
                    lambda i, v: logger.debug(f"{unit} n°{i}: {v}"),
                )

        return tqdm(
            iterable,
            desc=desc,
            leave=leave,
            disable=disable,
            unit=unit,
            **kwargs
        )
    

class GeneratorStep(Step, Generic[T]):
    """ A step that generate a dataset """
    output_dataset: Dataset[T] # Keeping generated dataset for logging purpose

    def __call__(self) -> Dataset[T]:
        from tqdm.contrib.logging import logging_redirect_tqdm
        with logging_redirect_tqdm():
            self.output_dataset = self.apply()
        return self.output_dataset
    
    @abstractmethod
    def apply(self) -> Dataset[T]:
        ...

    def __repr__(self) -> str:
        return f"{type(self).__name__}(output={self.output_dataset})"


class FilterStep(Step, Generic[T]):
    """ A step that filter a given dataset
    
    Thus, the __call__ method get a dataset in input and output a new one
    with the same kind of objects.

    The validated objects from the input dataset should always stay in the
    validated objects of the output dataset. Pending objects from the input
    dataset may be definitively rejected (not included in the output dataset),
    or remain in a pending state, or be definitively validated.
    """
    input_dataset: Dataset[T] # Keeping input dataset for logging purpose
    output_dataset: Dataset[T] # Keeping filtered dataset for logging purpose

    def __call__(self, dataset: Dataset[T], /) -> Dataset[T]:
        from tqdm.contrib.logging import logging_redirect_tqdm
        with logging_redirect_tqdm():
            self.input_dataset = dataset
            self.output_dataset = self.apply(self.input_dataset)
        return self.output_dataset
    
    @abstractmethod
    def apply(self, dataset: Dataset[T], /) -> Dataset[T]:
        ...

    def __repr__(self) -> str:
        return f"{type(self).__name__}(input={self.input_dataset}, output={self.output_dataset})"


class TransformerStep(Step, Generic[T, U]):
    """ A step that transforms a dataset into a dataset of different type """
    input_dataset: Dataset[T] # Keeping input dataset for logging purpose
    output_dataset: Dataset[U] # Keeping transformed dataset for logging purpose

    def __call__(self, dataset: Dataset[T], /) -> Dataset[U]:
        from tqdm.contrib.logging import logging_redirect_tqdm
        with logging_redirect_tqdm():
            self.input_dataset = dataset
            self.output_dataset = self.apply(self.input_dataset)
        return self.output_dataset

    @abstractmethod
    def apply(self, dataset: Dataset[T], /) -> Dataset[U]:
        ...

    def __repr__(self) -> str:
        return f"{type(self).__name__}(input={self.input_dataset}, output={self.output_dataset})"


###############################################################################
class GeneralStabilizerDimensionCheck(Step):
    """ Checking the dimension of the general stabilizer

    This step doesn't manage any dataset. It is just need to check that this
    package can be used with the given representation.

    It can actually be disabled using the no_dim_check option.
    """
    no_dim_check: bool

    def __init__(self, V: Representation, no_dim_check: bool = False, **kwargs: Any):
        super().__init__(V, **kwargs)
        self.no_dim_check = no_dim_check

    def apply(self) -> None:
        if self.no_dim_check:
            return
        
        from .stabK import dim_gen_stab_of_K
        Ms = self.V.actionK
        #MsR = [mat_C_to_R(M) for M in Ms.values()]
        # Check that the dim is computed in U_n(C)^s without the isolated S^1
        if (dim := dim_gen_stab_of_K(Ms)) > self.G.rank - self.V.dim_cone:
            raise ValueError(
                f"The general stabilizer of K in V is too big."
                f"Namely of dimension {dim}."
                f"The program does not work in this case."
            )

    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        group = parent_parser.add_argument_group(
            "Checking the dimension of the general stabilizer "
        )
        group.add_argument(
            "--no_dim_check",
            action="store_true",
            help="Don't check the dimension of the general stabilizer",
        )

    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "GeneralStabilizerDimensionCheck":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V=V,
            config=config,
            no_dim_check=config.no_dim_check,
            **kwargs,
        )


###############################################################################
class TauCandidatesStep(GeneratorStep[Tau]):
    """ A first list of dominent 1-PS whose kernel is supported at hyperplanes of weights
    
    It generates only pending Taus.
    """
    unique_tau: UniqueTauStr

    def __init__(
            self,
            V: Representation,
            unique_tau: UniqueTauStr = "SetOfTauCpp",
            **kwargs: Any
            ):
        super().__init__(V, **kwargs)
        # TODO when merged with dev_parallel2: compute so that 2^L > max_workers * chunk_size
        self.unique_tau = unique_tau

    def apply(self) -> Dataset[Tau]:
        from .tau import find_1PS
        return self.TDataset.from_separate(
            pending=self._tqdm(
                find_1PS(
                    self.V,
                    unique_tau=self.unique_tau,
                    quiet=self.quiet
                ),
                unit="tau"
            ),
            validated=[]
        )

    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        from typing import get_args

        group = parent_parser.add_argument_group(
            "First list of dominent 1-PS whose kernel is supported at hyperplanes of weights"
        )
        group.add_argument(
            "--unique_tau",
            type=lambda s: to_literal(UniqueTauStr, s),
            choices=get_args(UniqueTauStr),
            default="SetOfTauCpp",
            help="Type of the storage of unique tau",
        )

    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "TauCandidatesStep":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V=V,
            config=config,
            unique_tau=config.unique_tau,
            **kwargs,
        )


###############################################################################
class SubModuleConditionStep(FilterStep[Tau]):
    """
    Checking submodule condition
    
    It only reject pending Taus and doesn't modified the validated ones.
    """
    def apply(self, tau_dataset: Dataset[Tau]) -> Dataset[Tau]:
        from .parallel import Parallel
        from .utils import PartialFunction
        executor = Parallel().executor
        pending_tau = executor.filter(
            PartialFunction(Tau.is_sub_module, self.V),
            self._tqdm(tau_dataset.pending(), unit="tau"),
            chunk_size=executor.chunk_size * 32,
        )
        return self.TDataset.from_separate(
            pending=pending_tau,
            validated=tau_dataset.validated(),
        )
    

###############################################################################
class StabilizerConditionStep(FilterStep[Tau]):
    """
    Stabilizer condition
    
    It only reject pending Taus and doesn't modified the validated ones.
    """
    @staticmethod
    def tau_filter(tau: Tau, V: Representation) -> bool:
        from .stabK import dim_gen_stab_of_K

        Ms = V.actionK
        if  tau.is_dom_reg :
            return True
        else: 
            ListK=[beta.index_in_all_of_K(V.G) for beta in tau.orthogonal_rootsB]+[beta.opposite.index_in_all_of_K(V.G) for beta in tau.orthogonal_rootsU]
            ListChi=[V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)]+[V.dim+V.index_of_weight(chi) for chi in tau.orthogonal_weights(V)]
            return dim_gen_stab_of_K(Ms,ListK,ListChi) == V.G.rank - V.dim_cone + 1
        
    def apply(self, tau_dataset: Dataset[Tau]) -> Dataset[Tau]:
        from .parallel import Parallel
        from .utils import PartialFunction

        executor = Parallel().executor
        pending_tau = executor.filter(
            PartialFunction(StabilizerConditionStep.tau_filter, self.V),
            self._tqdm(tau_dataset.pending(), unit="tau"),
        )
        
        return self.TDataset.from_separate(
            pending=pending_tau,
            validated=tau_dataset.validated(),
        )


###############################################################################
class InequalityCandidatesStep(TransformerStep[Tau, Inequality]):
    """
    For each tau, computation the w with compatible tau-modules
    
    It generates only pending inequalities.
    """
    @staticmethod
    def List_Inv_Ws_Mod(tau: Tau, V: Representation) -> tuple[Tau, list[dict[int,list[Root]]]]:
        """ Helper method to avoid iterating two times on the tau """
        from .list_of_W import List_Inv_Ws_Mod
        return tau, List_Inv_Ws_Mod(tau, V)

    def apply(self, tau_dataset: Dataset[Tau]) -> Dataset[Inequality]:
        from .list_of_W import List_Inv_Ws_Mod
        from .parallel import Parallel
        from .utils import PartialFunction

        executor = Parallel().executor
        pending_tau = self._tqdm(tau_dataset.pending(), unit="tau")

        def ineq_generator() -> Iterator[Inequality]:
            inversions = executor.map(
                PartialFunction(InequalityCandidatesStep.List_Inv_Ws_Mod, self.V),
                pending_tau,
                chunk_size=executor.chunk_size * 2,
            )
            for tau, Lw in inversions:
                yield from (Inequality(tau, gr_inversions=gr_inv) for gr_inv in Lw)

        return self.TDataset.from_separate(
            pending=ineq_generator(),
            validated=[]
        )
    

###############################################################################
class TPiPreComputationStep(FilterStep[Inequality]):
    """
    Pre-computing the Representation.TPi3D matrix used in PiDominancy
    and Birationality steps.

    It thus filter nothing.
    """
    def apply(self, dataset: Dataset[T]) -> Dataset[T]:
        self.V.T_Pi_3D
        return dataset
    

###############################################################################
class PiDominancyStep(FilterStep[Inequality]):
    """
    Checking dominancy of the map pi
    
    It only reject pending inequalities and doesn't modified the validated ones.

    With the probabilistic method, valid inequalities may be rejected (with a low probability).
    """
    tpi_method: Method

    def __init__(self, V: Representation, tpi_method: Method = "probabilistic", **kwargs: Any):
        super().__init__(V, **kwargs)
        self.tpi_method = tpi_method

    def apply(self, ineq_dataset: Dataset[Inequality]) -> Dataset[Inequality]:
        from .list_of_W import Check_Rank_Tpi
        from .parallel import Parallel
        from .utils import PartialFunction

        executor = Parallel().executor
        inequalities = executor.filter(
            PartialFunction(Check_Rank_Tpi, self.V, self.tpi_method),
            self._tqdm(ineq_dataset.pending(), unit="ineq"),
            chunk_size=executor.chunk_size * 32,
        )
        return self.TDataset.from_separate(
            pending=inequalities,
            validated=ineq_dataset.validated(),
        )
    
    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        group = parent_parser.add_argument_group(
            "Checking dominancy of the map pi"
        )
        group.add_argument(
            "--tpi_method",
            type=lambda s: to_literal(Method, s),
            choices=("probabilistic", "symbolic"),
            default="probabilistic",
            help="Method for the surjectivity of Tpi",
        )

    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "PiDominancyStep":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V=V,
            config=config,
            tpi_method=config.tpi_method,
            **kwargs,
        )
    

###############################################################################
class LinearTriangularStep(FilterStep[Inequality]):
    """
    Checking Linear Triangular inequalities
    
    This filter can only definitively validate some of the inequalities (this inequalities are then not redondant).
    """
    def apply(self, ineq_dataset: Dataset[Inequality]) -> Dataset[Inequality]:
        from .linear_triangular import is_linear_triangular
        from itertools import chain
        pending: list[Inequality] = []
        validated: list[Inequality] = []

        pending_ineq = self._tqdm(ineq_dataset.pending(), unit="tau")
        def ineq_splitter() -> Iterator[tuple[Inequality, bool]]:
            for ineq in pending_ineq:
                if is_linear_triangular(self.V, ineq.tau, list(ineq.inversions)):
                    yield ineq, True # validated
                else:
                    yield ineq, False # still pending
        
        return self.TDataset.from_all(
            chain(
                ineq_splitter(),
                map(lambda ineq: (ineq, True), ineq_dataset.validated())
            )
        )

    

###############################################################################
class BKRConditionStep(FilterStep[Inequality]):
    """
    BKR condition
    
    It only reject pending inequalities and doesn't modified the validated ones.
    """
    kronecker: KroneckerCoefficient
    plethysm: PlethysmCache

    def __init__(self,
                 V: Representation,
                 kronecker: KroneckerCoefficient = KroneckerCoefficientMLCache(),
                 plethysm: PlethysmCache = PlethysmCache(),
                 **kwargs: Any):
        super().__init__(V, **kwargs)
        self.kronecker = kronecker
        self.plethysm = plethysm

    def apply(self, ineq_dataset: Dataset[Inequality]) -> Dataset[Inequality]:
        from .representation import ParticleRepresentation
        if isinstance(self.V, ParticleRepresentation) and self.G[0] >= 8:
            return self.TDataset.from_all(ineq_dataset.all())
        
        from .bkr import Multiplicity_SV_tau
        inequalities: list[Inequality] = []
        def ineq_filter(ineq: Inequality) -> bool:
            if list(ineq.inversions) == []:
                return True
            else:
                keep = Multiplicity_SV_tau(
                    ineq.tau,
                    ineq.weight_det(self.V),
                    self.V,
                    True,
                    self.kronecker, self.plethysm)
                assert isinstance(keep, bool)
                return keep
        
        return self.TDataset.from_separate(
            pending=filter(ineq_filter, self._tqdm(ineq_dataset.pending(), unit="ineq")),
            validated=ineq_dataset.validated(),
        )
    
    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        from typing import get_args
        method_choices = get_args(Method)

        group = parent_parser.add_argument_group(
            "BKR condition"
        )
        group.add_argument(
            "--kronecker",
            choices=("no_cache", "cache", "ml_cache"),
            default="ml_cache",
            help="Method for computing Kronecker coefficient"
        )

    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "BKRConditionStep":
        """ Build a step from the representation and the command-line arguments """
        from . import kronecker as kro

        kronecker: KroneckerCoefficient
        match config.kronecker:
            case "no_cache":
                kronecker = kro.KroneckerCoefficient()
            case "cache":
                kronecker = kro.KroneckerCoefficientCache()
            case "ml_cache":
                kronecker = kro.KroneckerCoefficientMLCache()
            case _:
                raise ValueError(f"Invalid Kronecker type {config.kronecker}")
            
        return super().from_config(
            V=V,
            config=config,
            kronecker=kronecker,
            plethysm=PlethysmCache(),
            **kwargs
        )
    

###############################################################################
class BirationalityStep(FilterStep[Inequality]):
    """
    Checking birationality (ramification divisor contracted) of the map pi
    
    It only outputs validated inequalities (non redondant) and reject all other inequalities
    (no pending inequalities remaining).
    """
    ram_schub_method: Method
    ram0_method: Method

    def __init__(self,
                 V: Representation,
                 ram_schub_method: Method = "probabilistic",
                 ram0_method: Method = "probabilistic",
                 **kwargs: Any):
        super().__init__(V, **kwargs)
        self.ram_schub_method = ram_schub_method
        self.ram0_method = ram0_method

    def apply(self, ineq_dataset: Dataset[Inequality]) -> Dataset[Inequality]:
        from .ramification import Is_Ram_contracted
        from .parallel import Parallel
        from .utils import PartialFunction
        from itertools import chain

        executor = Parallel().executor
        inequalities = executor.filter(
            PartialFunction(Is_Ram_contracted, self.V, self.ram_schub_method, self.ram0_method),
            self._tqdm(ineq_dataset.pending(), unit="ineq"),
        )

        return self.TDataset.from_separate(
            pending=[],
            validated=chain(inequalities, ineq_dataset.validated()),
        )
    
    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        from typing import get_args
        method_choices = get_args(Method)

        group = parent_parser.add_argument_group(
            "Checking birationality (ramification divisor contracted) of the map pi"
        )
        group.add_argument(
            "--ram_schub_method",
            type=lambda s: to_literal(Method, s),
            choices=method_choices,
            default="probabilistic",
            help="Method for checking if the Bruhat ramification divisors are contracted (p for probabilistic, s for symbolic)"
        )
        group.add_argument(
            "--ram0_method",
            type=lambda s: to_literal(Method, s),
            choices=method_choices,
            default="probabilistic",
            help="Method for checking if R0 is contracted (p for probabilistic, s for symbolic)"
        )
        
    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "BirationalityStep":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V=V,
            config=config,
            ram_schub_method=config.ram_schub_method,
            ram0_method=config.ram0_method,
            **kwargs,
        )
    

###############################################################################
class GrobnerStep(FilterStep[Inequality]):
    """
    Checking birationality via Grobner
        
    It only outputs validated inequalities (non redondant) and reject all other inequalities
    unless timeout is reached. In that case, the inequalities whose validation takes too much time
    are set as pending.
    """
    method: Method
    timeout: float

    def __init__(self,
                 V: Representation,
                 grobner_method: Method = "probabilistic",
                 grobner_timeout: float = 1,
                 **kwargs: Any):
        super().__init__(V, **kwargs)
        self.method = grobner_method
        self.timeout = grobner_timeout

    def apply(self, ineq_dataset: Dataset[Inequality]) -> ListDataset[Inequality]:
        from .groebner import Grobner_List_Test
        from itertools import chain
        grobner_true, grobner_inconclusive = Grobner_List_Test(
            self._tqdm(ineq_dataset.pending(), unit="ineq"),
            lim=self.timeout,
            V=self.V,
            method=self.method,
        )
        return ListDataset.from_separate(
            pending=grobner_inconclusive,
            validated=chain(ineq_dataset.validated(), grobner_true),
        )

    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        from typing import get_args
        method_choices = get_args(Method)

        group = parent_parser.add_argument_group(
            "Checking birationality via Grobner"
        )
        group.add_argument(
            "--grobner_method",
            type=lambda s: to_literal(Method, s),
            choices=method_choices,
            default="probabilistic",
            help="Method for checking if the fiber is a single point or not (p for probabilistic, s for symbolic)"
        )
        group.add_argument(
            "--grobner_timeout",
            type=float,
            default=1,
            help="Maximal processing time per inequality when checking birationaly (<= 0 to disable the limit)",
        )
        
    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "GrobnerStep":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V=V,
            config=config,
            grobner_method=config.grobner_method,
            grobner_timeout=config.grobner_timeout,
            **kwargs,
        )


###############################################################################
class ExportStep(FilterStep[Inequality]):
    """
    Exporting inequalities
    
    It exports all remaining inequalities (pending and validated).
    """
    formats: list[ExportFormat]

    def __init__(self,
                 V: Representation,
                 formats: ExportFormat | Iterable[ExportFormat] = [],
                 **kwargs: Any):
        super().__init__(V, **kwargs)
        if isinstance(formats, str):
            self.formats = [formats]
        else:
            self.formats = list(formats)

    def apply(self, ineq_dataset: Dataset[Inequality]) -> ListDataset[Inequality]:
        from .export import export_many
        inequations = ListDataset.from_all(
            self._tqdm(ineq_dataset.all(), unit="ineq")
        )
        export_many(self.formats, self.V, list(inequations))
        return inequations

    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        from typing import get_args
        method_choices = get_args(Method)

        group = parent_parser.add_argument_group(
            "Exporting inequalities"
        )
        group.add_argument(
            "--formats",
            type=lambda s: to_literal(ExportFormat, s),
            choices=get_args(ExportFormat),
            nargs='*',
            default=["Terminal"],
            help="Export format",
        )

    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "ExportStep":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V=V,
            config=config,
            formats=config.formats,
            **kwargs
        )


###############################################################################
InequalityFilterStr = Literal[
    "PiDominancy",
    "LinearTriangular",
    "BKRCondition",
    "Birationality",   
    "Grobner",
]

inequalities_filter_dict: Final[dict[InequalityFilterStr, type[Step]]] = {
    "PiDominancy": PiDominancyStep,
    "LinearTriangular": LinearTriangularStep,
    "BKRCondition": BKRConditionStep,
    "Birationality": BirationalityStep,
    "Grobner": GrobnerStep,
}

# Default filters
default_inequalities_filters: tuple[InequalityFilterStr, ...] = (
    "PiDominancy",
    "Birationality",
)

TStep = TypeVar("TStep", bound=Step)

class MomentConeStep(GeneratorStep[Inequality]):
    """
    Main step to generate inequalities of the moment cone
    
    It returns a dataset containing the inequalities that are definitively
    validated and the ones whose state is still pending.
    """
    config: Optional[Namespace] # Configuration from the command-line
    options: dict[str, Any] # Additional options passed to the constructor
    lazy: bool # Compute lazilly the inequalities without storing intermediate results
    store_steps: bool # Store all steps and their input/output datasets
    filters: list[InequalityFilterStr] # List of filters applied to the inequalities
    steps: list[Step] # All executed steps (for logging purpose)

    def __init__(
        self,
        V: Representation,
        filters: Iterable[str | InequalityFilterStr] = default_inequalities_filters,
        config: Optional[Namespace] = None,
        lazy: bool = False,
        store_steps: bool = False,
        **kwargs: Any,
    ):
        super().__init__(V, **kwargs)
        self.filters = [
            cast(InequalityFilterStr, to_literal(InequalityFilterStr, name))
            for name in filters
        ]
        self.config = config
        self.lazy = lazy
        self.store_steps = store_steps
        self.options = kwargs
        self.steps = []

    def __add_step(self, step_type: type[TStep]) -> TStep:
        """ Create and configure a new step """
        dataset_type = LazyDataset if self.lazy else ListDataset
        if self.config is None:
            step = step_type(self.V, dataset_type=dataset_type, **self.options)
        else:
            step = step_type.from_config(self.V, self.config, dataset_type=dataset_type)
        
        if self.store_steps:
            self.steps.append(step)

        return step
    
    def clear_steps(self) -> None:
        """ Clear all stored steps """
        self.steps.clear()
        
    def apply(self) -> Dataset[Inequality]:
        from .task import Task
        Task.quiet = self.quiet

        # Clearing TPi_3D cache to ensure using fresh random numbers
        self.V.clear_T_Pi_3D()

        # Clearing previous executed steps
        self.clear_steps()

        with Task(self.name) as main_task:
            # Checking if the cone has the expected dimension
            general_stab_dim_step = self.__add_step(GeneralStabilizerDimensionCheck)
            with Task(general_stab_dim_step.name):
                general_stab_dim_step()

            # Generate the list of candidates for tau
            tau_candidates: Dataset[Tau]
            tau_candidates_step = self.__add_step(TauCandidatesStep)
            with Task(tau_candidates_step.name) as task:
                tau_candidates = tau_candidates_step()
                task.log(f"tau_candidates: {tau_candidates}", indent=1)


            # Filters candidate tau
            for tau_filter_type in SubModuleConditionStep, StabilizerConditionStep:
                tau_filter_step = self.__add_step(tau_filter_type)
                with Task(tau_filter_step.name) as task:
                    tau_candidates = tau_filter_step(tau_candidates)
                    task.log(f"tau_candidates: {tau_candidates}", indent=1)
            
            # Transform tau to inequality
            ineq_candidates: Dataset[Inequality]
            ineq_candidates_step = self.__add_step(InequalityCandidatesStep)
            with Task(ineq_candidates_step.name) as task:
                ineq_candidates = ineq_candidates_step(tau_candidates)
                task.log(f"ineq_candidates: {ineq_candidates}", indent=1)

            # Pre-computation of Representation.TPi 3D matrix if necessary
            if 'PiDominancy' in self.filters or 'Birationality' in self.filters:
                TPi_step = self.__add_step(TPiPreComputationStep)
                with Task(TPi_step.name):
                    ineq_candidates = TPi_step(ineq_candidates)

            # Filters candidate inequalities
            for name in self.filters:
                ineq_filter_type = inequalities_filter_dict[name]
                ineq_filter_step = self.__add_step(ineq_filter_type)
                with Task(ineq_filter_step.name) as task:
                    ineq_candidates = ineq_filter_step(ineq_candidates)
                    task.log(f"ineq_candidates: {ineq_candidates}", indent=1)
            
            # Exporting inequalities
            export_step = self.__add_step(ExportStep)
            with Task(export_step.name):
                ineq_candidates = export_step(ineq_candidates)

            # Resume
            main_task.log("Steps resume:", indent=1)
            for step in self.steps:
                if isinstance(step, (GeneratorStep, FilterStep, TransformerStep)):
                    main_task.log(f"{step.name}: {step.output_dataset}", indent=2)

        
        return ineq_candidates

    @staticmethod
    def add_arguments(parent_parser: ArgumentParser, defaults: Mapping[str, Any] = {}) -> None:
        """ Add command-line arguments specific to this step """
        
        group = parent_parser.add_argument_group(
            "General moment cone computation options"
        )
        group.add_argument(
            "--quiet",
            action="store_true",
            help="Disable informations during computation",
        )
        group.add_argument( # FIXME: Move in step
            "--filters",
            type=lambda s: to_literal(InequalityFilterStr, s),
            nargs='*',
            choices=typing.get_args(InequalityFilterStr),
            default=default_inequalities_filters,
            help="Sequence of filters applied to the inequalities",
        )
        group.add_argument(
            "--lazy",
            action="store_true",
            help="Compute lazilly the inequalities (without storing intermediate results). When using --parallel, this option should be set together with --unordered.",
        )
        group.add_argument(
            "--store_steps",
            action="store_true",
            help="Store all intermediate steps and their input/ouptut datasets"
        )

        # Adding command-line options from other steps
        import sys
        module = sys.modules[__name__]
        for name in __all__:
            item_type = getattr(module, name)
            if issubclass(item_type, Step) and not issubclass(item_type, MomentConeStep):
                item_type.add_arguments(parent_parser)


    @classmethod
    def from_config(cls: type[Self], V: Representation, config: Namespace, **kwargs: Any) -> "MomentConeStep":
        """ Build a step from the representation and the command-line arguments """
        return super().from_config(
            V,
            config=config,
            filters=config.filters,
            lazy=config.lazy,
            store_steps=config.store_steps,
            **kwargs
        )
        

