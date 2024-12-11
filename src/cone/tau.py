from .typing import *
from .dimension import Dimension
from .blocks import Blocks
from .weight import Weight
from .root import Root
from .hyperplane_candidates import *
from .rings import matrix, ZZ, QQ
from .utils import extend_with_repetitions,Embeddings_mod_sym

import itertools
from functools import cached_property

"""
TODO:
- same method for both scalar product using @singledispatch?
- property or method for Pzero?
"""

__all__ = (
    "Tau",
    "ReducedTau",
    "find_1PS_reg_mod_sym_dim",
    "find_1PS_mod_sym_dim",
)

class Tau:
    """
    Tuple of partition along with a coefficient
    
    Example:
    >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)), 1)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1

    >>> d = Dimension((4, 3, 3))
    >>> tau = Tau.from_flatten((1, 3, 3, 2, 2, 2, 2, 1, 2, 2, 1), d)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1
    """
    #__slots__ = 'ccomponent', '_components' # FIXME cached_property cannot work without __dict__ ... => self managed cache or removing __slots__
    ccomponent: Optional[int]
    _components: Blocks[int]

    def __init__(self, components: Iterable[Sequence[int]] | Blocks[int], ccomponent: Optional[int] = None):
        """
        Tau initialization from a sequence of sub-group or directly from a partial matrix
        
        Beware that the ccomponent is at last position in the arguments since it is optional.
        """
        self.ccomponent = ccomponent
        if isinstance(components, Blocks):
            self._components = components.freeze()
        else:
            self._components = Blocks.from_blocks(components)

    @staticmethod
    def from_flatten(s: Iterable[int], d: Dimension) -> "Tau":
        """ Returns tau from its flattened version """
        all_components = tuple(s)
        if len(all_components) == d.sum:
            ccomponent = None
        elif len(all_components) == d.sum + 1:
            ccomponent = all_components[0]
            all_components = all_components[1:]
        else:
            raise ValueError("Invalid number of components")

        # Tau will be always immutable
        return Tau(Blocks.from_flatten(tuple(all_components), d), ccomponent)

    @staticmethod
    def from_zero_weights(S: Sequence[Weight], d: Dimension) -> "Tau":
        """ From a set of weights generating an hyperplane in X^*(T), returns a primitive Tau orthogonal to the hyperplane"""
        M = hyperplane_matrix(S, d)
        Mp=M.augment(matrix(ZZ,[len(d)*[0] for i in range(sum(d)+1)])) 
        for u in range(len(d)):
           shift=sum(d[k] for k in range(u))
           for i in range(d[u]):
               Mp[shift+i+1,u+len(S)]=1
        b=Mp.kernel().basis()
        if len(b)!=1:
           raise ValueError("Given set of weights does not generates an hyperplane")
        else:
           return Tau.from_flatten(b[0],d)

    def opposite(self) -> "Tau":
        """ Returns the opposite of some tau (same orthogonal hyperplane)
        """
        flatten_opp=[-self.ccomponent]
        dd=[]
        for comp in self.components:
           flatten_opp+=[-x for x in comp]
           dd.append(len(comp))
        return Tau.from_flatten(flatten_opp,Dimension(dd))

    def __len__(self) -> int:
        """ Number of components """
        return len(self._components)

    @cached_property
    def d(self) -> Dimension:
        """ Length of each component """
        return Dimension(self._components.sizes)

    @cached_property
    def components(self) -> tuple[tuple[int, ...], ...]:
        """ Sequence of the components of tau """
        return cast(tuple[tuple[int, ...], ...], tuple(self._components.blocks))

    @property
    def flattened(self) -> Iterable[int]:
        """
        Returns the whole tau as a unique sequence.

        ie (cc-component | column1 | column2 | ...)
        """
        if self.ccomponent is None:
            return self._components.flatten
        else:
            return (self.ccomponent,) + cast(tuple[int], self._components.flatten)

    @cached_property
    def reduced(self) -> "ReducedTau":
        """ Returns reduced form of tau """
        return ReducedTau(self)
        
    def m_extend_with_repetitions(self,d:Dimension)-> Iterable["Tau"]:
        extend_each_comp=[extend_with_repetitions(x,d[i]) for i,x in enumerate(self._components)]
        for ext in itertools.product(*extend_each_comp):
            yield Tau(ext,self.ccomponent)

    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + repr(self._components)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Tau):
            return NotImplemented
        return self.ccomponent == other.ccomponent and self._components == other._components

    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash((self.ccomponent, self._components))

    # TODO: or same dot method with @singledispatch ?
    def dot_weight(self, weight: Weight, ccomponent: Optional[int] = None) -> int:
        """ Scalar product of tau with a weight of V """
        assert not (self.ccomponent is None and ccomponent is None)
        if ccomponent is None:
            ccomponent = self.ccomponent
        return cast(int, ccomponent) + sum(c[wi] for c, wi in zip(self.components, weight))

    def dot_root(self, root: Root) -> int:
        """ Scalar product of tau with a root of U """
        c = self.components[root.k]
        return c[root.i] - c[root.j]

    @cached_property
    def is_dom_reg(self) -> bool:
        """ Check if tau is dominant and regular """
        return all(all(a > b for a, b in itertools.pairwise(c)) for c in self.components)

    @cached_property
    def is_dominant(self) -> bool:
        """ Check if tau is dominant """
        return all(all(a >= b for a, b in itertools.pairwise(c)) for c in self.components)

    @cached_property
    def sl_representative(self) -> "Tau":
        """ Returns representative of tau in C^* x (SLn)^3 """
        assert self.ccomponent is not None
        from math import lcm, gcd
        tau_lcm = lcm(*self.flattened)
        ccomponent = self.ccomponent * tau_lcm
        columns = []
        for dj, cj in zip(self.d, self.components):
            column_sum = sum(cj)
            shift = column_sum * tau_lcm // dj
            columns.append([tau_lcm * cji - shift for cji in cj])
            ccomponent += shift

        res_gcd = gcd(ccomponent, *itertools.chain.from_iterable(columns))
        return Tau(
            tuple(tuple(v // res_gcd for v in cj) for cj in columns),
            ccomponent // res_gcd
        )

    @cached_property
    def end0_representative(self) -> "Tau":
        """ Returns representative of tau in ??? space (final value of each block is zero) """
        assert self.ccomponent is not None
        from math import gcd
        total_shift = 0
        columns = []
        for cj in self.components:
            columns.append([cji - cj[-1] for cji in cj])
            total_shift += cj[-1]
        ccomponent = self.ccomponent + total_shift

        res_gcd = gcd(ccomponent, *itertools.chain.from_iterable(columns))
        return Tau(
            tuple(tuple(v // res_gcd for v in cj) for cj in columns),
            ccomponent // res_gcd
        )

    def grading_dictionary(self, basis:list[T], dot_basis) -> dict[int,list[T]]:
        """ basis is a set of elements parameterizing a basis of the space of interest made of eigenvector under action of tau.
        dot_basis T -> int is a function associating to an element of the basis its eigenvalue.
        grading_dictionary returns a dictionary whose keys are eigenvalues. For each key p, the values in the entry p correspond to a basis of the eigenspace
        """
        result: dict[int, list[T]] = {}
        for chi in basis:
            p = dot_basis(chi)
            result.setdefault(p, []).append(chi)
        return result

    @cached_property
    def grading_weights(self) -> dict[int, list[Weight]]:
        """ dictionary whose keys are eigenvalues of the action of tau on V. For each key p, the weights in the entry p correspond to a basis of the eigenspace """
        weights = Weight.all(self.d)
        return self.grading_dictionary(weights,self.dot_weight)
        # Former version
        #weights = Weight.all(self.d)
        #result: dict[int, list[Weight]] = {}
        #for chi in weights:
        #    p = self.dot_weight(chi)
        #    result.setdefault(p, []).append(chi)
        #return result
        

    @cached_property
    def grading_roots(self) -> dict[int, list[Root]]:
        """ dictionary whose keys are eigenvalues of the action of tau on u (sum of positive root spaces). For each key p, the roots in the entry p correspond to a basis of the eigenspace """
        roots = Root.all_of_U(self.d)
        return self.grading_dictionary(roots,self.dot_root)
    
    def filter_dict(self, dic, prop) -> dict[int, list]:
        """ Selects in the dictionary dic, the keys satisfying the property prop"""
        def property_key(pair):
           x,v=pair
           return prop(x)
        return dict(filter(property_key, dic.items()))

    @property
    def positive_weights(self) -> dict[int, list[Weight]]:
        return self.filter_dict(self.grading_weights,lambda x: x>0)
    
    @property
    def non_negative_weights(self) -> dict[int, list[Weight]]:
        return self.filter_dict(self.grading_weights,lambda x: x>=0)
    
    @property
    def positive_roots(self) -> dict[int, list[Root]]:
        return self.filter_dict(self.grading_roots,lambda x: x>0)

    #TODO: definition changed thanks to above definition.
    #Optimization might be done to avoid lists and dict, dealing only with iterables?
    @property
    def orthogonal_roots(self) -> Iterable[Root]:
        """ All the root beta so that <beta, tau> = 0 """
        if 0 in self.grading_roots:
           return self.grading_roots[0]
        else: 
           return []

    @property
    def orthogonal_weights(self) -> Iterable[Root]:
        """ All the root beta so that <beta, tau> = 0 """
        if 0 in self.grading_weights:
           return self.grading_weights[0]
        else: 
           return []

    @cached_property
    def sort_mod_sym_dim(self) -> "Tau":
        """
        Sort tau by block of the dimensions
        
        >>> d = Dimension((2, 2, 2, 1, 1, 1))
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], d)
        >>> tau
        1 | 6 2 | 1 4 | 1 4 | 5 | 3 | 1
        >>> tau.sort_mod_sym_dim
        1 | 1 4 | 1 4 | 6 2 | 1 | 3 | 5
        """
        blocks = (sorted(b) for b in Blocks(self.components, self.d.symmetries))
        return Tau(itertools.chain.from_iterable(blocks), self.ccomponent)
        
    def orbit_symmetries(self) -> Iterable["Tau"]:
        """
        Lists the orbit of tau under symmetries of dimensions of its components  
        """
        for sym_comp in self._components.orbit_symmetries():
           yield Tau(sym_comp,self.ccomponent)

    @cached_property
    def dim_Pu(self) -> int:
        """
        Dimension of Pu
        
        Example:
        >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)), 1)
        >>> tau
        1 | 3 3 2 2 | 2 2 1 | 2 2 1
        >>> tau.dim_Pu
        8
        >>> sum(len(roots) for roots in tau.positive_roots.values())
        8
        """
        sum_di2 = sum(di**2 for di in self.d)
        sum_mi2 = sum(mi**2 for mi in self.reduced.mult.flatten)
        return (sum_di2 - sum_mi2) // 2

class ReducedTau:
    """
    Tau in a reduction form
    
    Example:
    >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)), 1)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1
    >>> rtau = tau.reduced
    >>> rtau
    1 | 3^2 2^2 | 2^2 1^1 | 2^2 1^1
    >>> rtau.values
    3 2 | 2 1 | 2 1
    >>> rtau.mult
    2 2 | 2 1 | 2 1
    """
    __slots__ = 'ccomponent', 'values', 'mult'
    ccomponent: Optional[int]
    values: Blocks[int]
    mult: Blocks[int]

    def __init__(self, tau: Tau):
        from .utils import group_by_block
        values = Blocks.from_flatten([], (0,) * len(tau.d))
        mult = Blocks.from_flatten([], (0,) * len(tau.d))

        for j, component in enumerate(tau.components):
            values[j], mult[j] = zip(*group_by_block(component))

        self.ccomponent = tau.ccomponent
        self.values = values.freeze()
        self.mult = mult.freeze()

    def __len__(self) -> int:
        """ Number of components """
        return len(self.values)

    @property
    def small_d(self) -> Dimension:
        return Dimension(self.values.sizes)

    # FIXME: it seems that i,j is not really needed => simplification?
    def __getitem__(self, idx: tuple[int, int]) -> tuple[int, int]:
        return self.values[idx[1], idx[0]], self.mult[idx[1], idx[0]]

    def __repr__(self) -> str:
        return f"{self.ccomponent} | " + " | ".join(
            " ".join(f"{v}^{m}" for v, m in zip(cv, cm))
            for cv, cm in zip(self.values.blocks, self.mult.blocks)
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ReducedTau):
            return NotImplemented
        return self.ccomponent == other.ccomponent and self.values == other.values and self.mult == other.mult

    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash((self.ccomponent, self.values, self.mult))

    # FIXME: returns an generator and not a sequence. Is it adapted to a property (or a cached_property?)
    # TODO: it will be used many time => returns a sequence and @cached_property
    # TODO: remove? information already contained in orthogonal weight like subsequent function.
    @property
    def Pzero(self) -> Iterable[Weight]:
        """ Search for weights w of V so that C_component + sum_k tau_red[w_k, k] = 0 """
        as_tau = Tau(self.values, self.ccomponent) # We use the scalar product of Tau
        for weight in Weight.all(self.small_d):
            if as_tau.dot_weight(weight) == 0:
                yield weight

def find_1PS_reg_mod_sym_dim(d:Dimension,u) -> Sequence["Tau"]:
    """
    Lists candidates dominant regular 1 parameter subgroups tau, whose orthogonal is generated by weights, and such that there is at most u weights we of V with tau(we)>0
    With results up to the action of the symmetries of d.
    """
    Liste_hr=find_hyperplanes_reg_mod_sym_dim(d,u)
    Liste_1PS=list(set([Tau.from_zero_weights(h,d).end0_representative.sort_mod_sym_dim for h in Liste_hr]))
    Liste_1PS_sign=Liste_1PS+[tau.opposite() for tau in Liste_1PS]
    return([tau for tau in Liste_1PS_sign if tau.is_dom_reg])

def find_1PS_mod_sym_dim(d: Dimension) -> Sequence["Tau"]:
    """
    Same as find_1PS_reg_mod_sym_dim without regularity condition
    Computed by 
    """
    # Initialisation with regular 1-PS
    Liste_1PS=find_1PS_reg_mod_sym_dim(d,d.uMAX(d))
    print('For d=',d,'we get',len(Liste_1PS),' candidates regular dominant up to symmetry')
    # Looking for 1-PS by extension
    sub_dim=[Dimension(p) for p in Partition(d).all_subpartitions()][1:-1] #[1:-1] excludes 1... 1 and d
    for small_d in sub_dim:  
        umax=uMAX(d,small_d)
        #Recover by induction all candidates 1-PS mod symmetry
        Liste_1PS_smalld_mod_sym= hyperplane_reg_mod_sym_dim(small_d,umax)
        print('For d=',small_d,'we get',len(Liste_1PS_small_mod_sym),' candidates regular dominant up to symmetry')
        Liste_1PS_smalld=sum([list(tau.orbit_symmetries()) for tau in Liste_1PS_smalld_mod_sym]  ,[])
        #from now on, the program has to be adapted
        for permut in Embeddings_mod_sym(d,small_d):
            for tau in Liste_1PS_smalld:
                tau_twist=Tau([tau0._components[i] for i in permut],tau0.ccomponent)
                list_tau_extended=tau_twist.m_extend_with_repetitions(d)
                #list_tau_extended=[[tau[0]]+l for l in list_tau_extended]
                # Check dimU
                list_tau_extended_dimU=[]
                for tau_ext in list_tau_extended:
                    if Istaunlequ(tau_ext,d,dimU_GL_power_s(tau_ext,d),Poids):
                        list_tau_extended_dimU.append(tau_ext)
                Liste+=list_tau_extended_dimU
    return(Liste)
