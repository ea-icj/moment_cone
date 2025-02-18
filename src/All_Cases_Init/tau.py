#from typing import Generic, Any

from .typing import *
from .utils import *
from .permutation import Permutation
from .partition import *
from .group import *
from .weight import *
from .rep import *
from .hyperplane_candidates import *
from .blocks import Blocks
from .root import Root
from .rings import matrix, ZZ, QQ

#from collections.abc import Iterable, Sequence, Callable


import itertools
from functools import cached_property


__all__ = (
    "Tau",
    "ReducedTau",
    "unique_modulo_symmetry_list_of_tau",
    "full_under_symmetry_list_of_tau",
    "find_1PS",
)

class Tau:
    """
    Tuple of partition along with a coefficient
    
    Example:
    >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)))
    >>> tau
    3 3 2 2 | 2 2 1 | 2 2 1
    >>> tau.G
    GL(4)xGL(3)xGL(3)
    >>> G = LinGroup([1, 4, 3, 3])
    >>> tau = Tau.from_flatten((1, 3, 3, 2, 2, 2, 2, 1, 2, 2, 1), G)
    >>> tau
    1 | 3 3 2 2 | 2 2 1 | 2 2 1
    """
    #__slots__ =  '_components' # FIXME cached_property cannot work without __dict__ ... => self managed cache or removing __slots__
    _components: Blocks[int]

    def __init__(self, components: Iterable[Sequence[int]] | Blocks[int]):
        """
        Tau initialization from a sequence of sub-group or directly from a partial matrix
        
        """
        
        if isinstance(components, Blocks):
            self._components = components.freeze()
        else:
            self._components = Blocks.from_blocks(components)
        #self.G = LinGroup([len(c) for c in components])
        
    @staticmethod
    def from_flatten(s: Iterable[int], G: LinGroup) -> "Tau":
        """ Returns tau from its flattened version """
        all_components = tuple(s)
        if len(all_components) != G.rank:
            raise ValueError("Invalid number of components")

        # Tau will be always immutable
        return Tau(Blocks.from_flatten(tuple(all_components), G))
    
    @cached_property
    def G(self) -> LinGroup:
        return LinGroup([len(c) for c in self.components])

    @staticmethod
    def from_zero_weights(weights: Sequence[Weight], V: Representation) -> "Tau":
        """
        From a set of weights generating an hyperplane in X^*(T), returns a primitive Tau orthogonal to the hyperplane

        TODO: doctest
        """

        M = matrix([v.as_vector.list() for v in weights]).transpose()
        if V.type=='kron': # In this case we add conditions of normalization because tau is defined in X_*(T/Z)
                           # We choose each block ends by 0
            M = M.augment(matrix(ZZ, [(len(V.G)-1) * [0] for i in range(V.G.rank)]))
            shift_i = 0
            shift_j = len(weights)
            for j,d in enumerate(V.G[:-1]):           
               M[shift_i + d-1,shift_j+j] = 1
               shift_i += d

        b = M.kernel().basis()
        if len(b) != 1:
           raise ValueError("Given set of weights does not generates an hyperplane")
        else:
           return Tau.from_flatten(b[0], V.G)

    @property
    def opposite(self) -> "Tau":
        """
        Returns the opposite of some tau (same orthogonal hyperplane)

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)))
        >>> tau
        3 2 2 | 4 2 1 | 3 2
        >>> tau.opposite
        -3 -2 -2 | -4 -2 -1 | -3 -2
        """
        return Tau.from_flatten(
            (-x for x in self.flattened),
            self.G
        )

    def __len__(self) -> int:
        """ Number of components """
        return len(self.G)

    @cached_property
    def components(self) -> tuple[tuple[int, ...], ...]:
        """ Sequence of the components of tau """
        return cast(tuple[tuple[int, ...], ...], tuple(self._components.blocks))

    @property
    def flattened(self) -> Iterable[int]:
        """
        Returns the whole tau as a unique sequence.

        ie ( column1 | column2 | ...)
        """
        return self._components.flatten
      
    @cached_property
    def reduced(self) -> "ReducedTau":
        """ 
        Returns reduced form of tau
        
        Example:
        >>> tau = Tau(((3, 3, 2), (4, 4, 2), (1, 1)))
        >>> tau.reduced
        3^2 2^1 | 4^2 2^1 | 1^2
        """
        return ReducedTau(self)
    
    def extend_from_S(self, p: Partition):
        my_list=[]
        for i,x in zip(p,self.flattened):
            my_list+=i*[x]
        my_list.sort(reverse=True)
        return Tau.from_flatten(my_list,LinGroup([sum(p)]))
    
    def m_extend_with_repetitions(self,G: LinGroup)-> Iterable["Tau"]:
        """
        Extends self in all possible manners to a Tau of G by adding repetitions
                
        Example:
        >>> tau = Tau([[5, 4, 0], [4, 3, 0], [4, 2, 0]])
        >>> tau
        5 4 0 | 4 3 0 | 4 2 0
        >>> for t in tau.m_extend_with_repetitions(LinGroup([4, 4, 3])):
        ...     print(t)
        5 4 0 0 | 4 3 0 0 | 4 2 0
        5 4 0 0 | 4 3 3 0 | 4 2 0
        5 4 0 0 | 4 4 3 0 | 4 2 0
        5 4 4 0 | 4 3 0 0 | 4 2 0
        5 4 4 0 | 4 3 3 0 | 4 2 0
        5 4 4 0 | 4 4 3 0 | 4 2 0
        5 5 4 0 | 4 3 0 0 | 4 2 0
        5 5 4 0 | 4 3 3 0 | 4 2 0
        5 5 4 0 | 4 4 3 0 | 4 2 0
        """
        extend_each_comp=[extend_with_repetitions(x,G[i]) for i,x in enumerate(self.components)]
        for ext in itertools.product(*extend_each_comp):
            yield Tau(ext)

    def __repr__(self) -> str:
        return repr(self._components)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, Tau):
            return NotImplemented
        return  self._components == other._components

    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash(self._components)

    # TODO: or same dot method with @singledispatch ?
    def dot_weight(self, chi ) -> int: # chi can be a Weight or a Weight_as_list or  a Weight_as_list_of_list 
        """ Scalar product of tau with a weight of V """
        if chi.as_list is not None : 
            return sum(c[eps] for c, eps in zip(self.components, chi.as_list))
        if chi.as_list_of_list is not None :
            tot=0
            for c,eps in zip(self.components,chi.as_list_of_list):
                tot+=sum([c[i] for i in eps])
            return tot
        return sum([x*y for x,y in zip(self.flattened,chi.as_vector)])
    
    def dot_root(self, root: Root) -> int:
        """ 
        Scalar product of tau with a root of self.G
        
        Example:
        >>> tau = Tau([[5, 4, 0], [4, 3, 0], [4, 2, 0]])
        >>> tau.dot_root(Root(1,0,2))
        4
        """
        c = self.components[root.k]
        return c[root.i] - c[root.j]

    @cached_property
    def is_dom_reg(self) -> bool:
        """ 
        Check if tau is dominant and regular 
        Examples:
        >>> tau=Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)))
        >>> tau.is_dom_reg
        False
        >>> tau=Tau(((3, 2, 1), (2,  1), (2, 1)))
        >>> tau.is_dom_reg
        True

        """
        return all(all(a > b for a, b in itertools.pairwise(c)) for c in self.components)

    @cached_property
    def is_dominant(self) -> bool:
        """ 
        Check if tau is dominant 
        
        Example:
        >>> tau=Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1)))
        >>> tau.is_dominant
        True
        """
        return all(all(a >= b for a, b in itertools.pairwise(c)) for c in self.components)
    

    def grading_weights_in(self, weights: Iterable[Weight]) -> dict[int, list[Weight]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on a given subset of V.
        """
        from .utils import grading_dictionary
        return grading_dictionary(weights, self.dot_weight)

    
    def grading_weights(self,V : Representation) -> dict[int, list[Weight]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on V.
        
        For each key p, the weights in the entry p correspond to a basis of the eigen space
        
        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2),(1,)))
        >>> tau
        3 2 2 | 4 2 1 | 3 2 | 1
        >>> V=Reprensentation(tau.G,'kron')
        >>> gw = tau.grading_weights(V)
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        -2: [Weight((1, 2, 1), idx: 11), Weight((2, 2, 1), idx: 17)]
        -1: [Weight((0, 2, 1), idx: 5), Weight((1, 1, 1), idx: 9), Weight((1, 2, 0), idx: 10), Weight((2, 1, 1), idx: 15), Weight((2, 2, 0), idx: 16)]
        0: [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]
        1: [Weight((0, 1, 0), idx: 2), Weight((1, 0, 1), idx: 7), Weight((2, 0, 1), idx: 13)]
        2: [Weight((0, 0, 1), idx: 1), Weight((1, 0, 0), idx: 6), Weight((2, 0, 0), idx: 12)]
        3: [Weight((0, 0, 0), idx: 0)]
        """
        return self.grading_weights_in(V.all_weights)

    def grading_roots_in(self, roots: Iterable[Root]) -> dict[int, list[Root]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on the given subset of U.
        """
        from .utils import grading_dictionary
        return grading_dictionary(roots, self.dot_root)


    @cached_property
    def grading_rootsU(self) -> dict[int, list[Root]]:
        """
        Dictionary whose keys are eigenvalues of the action of tau on U (sum of positive root spaces).
        
        For each key p, the roots in the entry p correspond to a basis of the eigen space.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gr = tau.grading_roots
        >>> for k in sorted(gr.keys()):
        ...     print(f"{k}:", gr[k])
        0: [Root(k=0, i=1, j=2)]
        1: [Root(k=0, i=0, j=1), Root(k=0, i=0, j=2), Root(k=1, i=1, j=2), Root(k=2, i=0, j=1)]
        2: [Root(k=1, i=0, j=1)]
        3: [Root(k=1, i=0, j=2)]
        """
        return self.grading_roots_in(Root.all_of_U(self.G))

    @cached_property
    def grading_rootsB(self) -> dict[int, list[Root]]:
        return self.grading_roots_in(Root.all_of_B(self.G))
        
    def positive_weights(self, V: Representation) -> dict[int, list[Weight]]:
        """
        Basis of the eigen space for positive eigen values for the action of tau on V.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.positive_weights
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        1: [Weight((0, 1, 0), idx: 2), Weight((1, 0, 1), idx: 7), Weight((2, 0, 1), idx: 13)]
        2: [Weight((0, 0, 1), idx: 1), Weight((1, 0, 0), idx: 6), Weight((2, 0, 0), idx: 12)]
        3: [Weight((0, 0, 0), idx: 0)]
        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_weights(V), lambda x: x > 0)

    
    def non_positive_weights(self,V: Representation) -> dict[int, list[Weight]]:
        """
        Basis of the eigen space for non-positive eigen values for the action of tau on V.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.non_positive_weights
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        -2: [Weight((1, 2, 1), idx: 11), Weight((2, 2, 1), idx: 17)]
        -1: [Weight((0, 2, 1), idx: 5), Weight((1, 1, 1), idx: 9), Weight((1, 2, 0), idx: 10), Weight((2, 1, 1), idx: 15), Weight((2, 2, 0), idx: 16)]
        0: [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]

        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_weights(V), lambda x: x <= 0)


    
    def non_negative_weights(self,V: Representation) -> dict[int, list[Weight]]:
        """
        Basis of the eigen space for non-negative eigen values for the action of tau on V.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> gw = tau.non_negative_weights(V)
        >>> for k in sorted(gw.keys()):
        ...     print(f"{k}:", gw[k])
        0: [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]
        1: [Weight((0, 1, 0), idx: 2), Weight((1, 0, 1), idx: 7), Weight((2, 0, 1), idx: 13)]
        2: [Weight((0, 0, 1), idx: 1), Weight((1, 0, 0), idx: 6), Weight((2, 0, 0), idx: 12)]
        3: [Weight((0, 0, 0), idx: 0)]
        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_weights(V), lambda x: x >= 0)
    
    @property
    def positive_rootsU(self) -> dict[int, list[Root]]:
        """
        Basis of the eigen space for positive eigen values for the action of tau on U.

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2),(-7,)))
        >>> tau
        3 2 2 | 4 2 1 | 3 2 | -7
        >>> gr = tau.positive_roots
        >>> for k in sorted(gr.keys()):
        ...     print(f"{k}:", gr[k])
        1: [Root(k=0, i=0, j=1), Root(k=0, i=0, j=2), Root(k=1, i=1, j=2), Root(k=2, i=0, j=1)]
        2: [Root(k=1, i=0, j=1)]
        3: [Root(k=1, i=0, j=2)]
        """
        from .utils import filter_dict_by_key
        return filter_dict_by_key(self.grading_rootsU, lambda x: x > 0)

    #TODO: definition changed thanks to above definition.
    #Optimization might be done to avoid lists and dict, dealing only with iterables?
    @property
    def orthogonal_rootsU(self) -> list[Root]:
        """
        All the positive roots beta of V so that <beta, tau> = 0

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> tau.orthogonal_roots
        [Root(k=0, i=1, j=2)]
        """
        return self.grading_rootsU.get(0, [])

    @property
    def orthogonal_rootsB(self) -> list[Root]:
        """
        All the positive of zero roots beta of V so that <beta, tau> = 0
        """
        return self.grading_rootsB.get(0, [])
    
    @cached_property
    def orthogonal_rootsK(self) -> list[Root]:
        """
        All the roots beta of V so that <beta, tau> = 0. beta=Root(k,i,i) allowed. 
        We get a set indexing a bases of Lie(K^tau)

        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> for t in tau.orthogonal_all_roots:
        ...     print(t)
        Root(k=0, i=0, j=0)
        Root(k=0, i=1, j=1)
        Root(k=0, i=1, j=2)
        Root(k=0, i=2, j=1)
        Root(k=0, i=2, j=2)
        Root(k=1, i=0, j=0)
        Root(k=1, i=1, j=1)
        Root(k=1, i=2, j=2)
        Root(k=2, i=0, j=0)
        Root(k=2, i=1, j=1)
        """
        res=[]
        for alpha in Root.all_of_K(self.G):
            if self.dot_root(alpha)==0:
                res.append(alpha)
        return(res)        

    
    def orthogonal_weights(self,V: Representation) -> list[Weight]:
        """
        All the root beta so that <beta, tau> = 0
        
        >>> tau = Tau(((3, 2, 2), (4, 2, 1), (3, 2)), -7)
        >>> tau
        -7 | 3 2 2 | 4 2 1 | 3 2
        >>> tau.orthogonal_weights(V)
        [Weight((0, 1, 1), idx: 3), Weight((0, 2, 0), idx: 4), Weight((1, 1, 0), idx: 8), Weight((2, 1, 0), idx: 14)]
        """
        return self.grading_weights(V).get(0, [])

    @cached_property
    def sort_mod_sym_dim(self) -> "Tau":
        """
        Sort tau by block of the dimensions
        
        >>> G = LinGroup([2, 2, 2, 1, 1, 1])
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3, 1], G)
        >>> tau
        1 | 6 2 | 1 4 | 1 4 | 5 | 3 | 1
        >>> tau.sort_mod_sym_dim
        1 | 1 4 | 1 4 | 6 2 | 1 | 3 | 5
        """
        blocks = (sorted(b) for b in Blocks(self.components, self.G.outer))
        return Tau(itertools.chain.from_iterable(blocks))
        
    def orbit_symmetries(self) -> Iterable["Tau"]:
        """
        Lists the orbit of tau under symmetries of dimensions of its components

        Example:
        >>> G = LinGroup([2, 2, 2, 1, 1])
        >>> tau = Tau.from_flatten([1, 6, 2, 1, 4, 1, 4, 5, 3], G)
        >>> tau
        1 | 6 2 | 1 4 | 1 4 | 5 | 3
        >>> for t in tau.orbit_symmetries():
        ...     print(t)
        1 | 1 4 | 1 4 | 6 2 | 3 | 5
        1 | 1 4 | 1 4 | 6 2 | 5 | 3
        1 | 1 4 | 6 2 | 1 4 | 3 | 5
        1 | 1 4 | 6 2 | 1 4 | 5 | 3
        1 | 6 2 | 1 4 | 1 4 | 3 | 5
        1 | 6 2 | 1 4 | 1 4 | 5 | 3
        """
        for sym_comp in orbit_symmetries(self._components, self.G.outer):
            yield Tau(sym_comp)

    #@staticmethod
    def is_sub_module(self,V: Representation) -> bool :
        """
        Check if V^{\tau>0} is a C^*-submodule of U

        TODO: doctest
        """
        # Create the dictionnary for V^{\tau>0}
        gw = self.positive_weights(V) 
        # Create the dictionnary for U
        # FIXME: it's self.grading_roots
        gr = grading_dictionary(Root.all_of_U(self.G), self.dot_root)

        # Compare them
        # TODO: using all instead (after doctest) ?
        for p in gw.keys():
            if p not in gr.keys() or len(gr[p])<len(gw[p]) :
                return False
        return True

    @cached_property
    def dim_Pu(self) -> int:
        """
        Dimension of Pu
        
        Example:
        >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1),(1,)))
        >>> tau
        3 3 2 2 | 2 2 1 | 2 2 1 | 1
        >>> tau.dim_Pu
        8
        >>> sum(len(roots) for roots in tau.positive_roots.values())
        8
        """
        sum_di2 = self.G.dim
        sum_mi2 = sum(mi**2 for mi in self.reduced.mult.flatten)
        return (sum_di2 - sum_mi2) // 2

    @cached_property
    def end0_representative(self) -> "Tau":
        """ 
        Returns representative of tau in X*(T/Z)  with final value of each block is zero.
        Since Z is trivial when len(G)==1 return self.

        Examples:
        >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1),(1,)))
        >>> tau.end0_representative
        1 1 0 0 | 1 1 0 | 1 1 0 | 5
        >>> tau = Tau(((3, 3, 2, 2),))
        >>> tau.end0_representative
        3 3 2 2
        """

        # Return self if not kron
        if len(self.G)==1:
            return self

        from sage.all import gcd
        total_shift = 0
        columns = []
        for cj in self.components[:-1]: # Components excepted that of the last GL(1)
            columns+=[cji - cj[-1] for cji in cj]
            total_shift += cj[-1]
        last_component = self.flattened[-1] + total_shift
        columns.append(last_component)
        res_gcd = gcd(columns)
        
        return Tau.from_flatten([x // res_gcd for x in columns], self.G)
    
    @cached_property
    def sl_representative(self) -> "Tau":
        """ 
        Returns representative of tau in X*(T/Z)  in product of SL.
        Since Z is trivial when len(G)==1 return self.

        Examples:
        >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1),(1,)))
        >>> tau.sl_representative
        3 3 -3 -3 | 2 2 -4 | 2 2 -4 | 41
        >>> tau = Tau(((3, 3, 2, 2),))
        >>> tau.sl_representative
        3 3 2 2
        """
        if len(self.G)==1:
            return self
        from math import lcm, gcd
        #flat_non_zero=[x for x in self.flattened if x!=0]
        d=(len(v) for v in self.components)
        tau_lcm = lcm(*d)
        ccomponent = self.components[-1][0] * tau_lcm
        columns = []
        for cj in self.components[:-1]:
            dj=len(cj)
            column_sum = sum(cj)
            shift = column_sum * tau_lcm // dj
            columns.append([tau_lcm * cji - shift for cji in cj])
            ccomponent += shift
        res_gcd = gcd(ccomponent, *itertools.chain.from_iterable(columns))
        return Tau(
            tuple(tuple(v // res_gcd for v in cj) for cj in columns) +((ccomponent // res_gcd,),)
        )


    def summands_Vtau(self,V : Representation)  -> list[list[int]] : 
        """
        V^tau can be written as a direct sum. Compute the index set of this sum. Namely :
        For Kron : Create the list of lists L=[i_0, i_2, ..., i_{s-1}] such that  sum_k tau_red[k][i_k] = 0
        For Fermion or Boson : s=len(tau_red.values). Create the list of lists L=[i_0, i_2, ..., i_{s-1}] such that  sum_k i_k * tau_red[k] = 0
        """
        #tau_red=Tau(self.reduced.values,LinGroup([len(x) for x in self.reduced.values]))
        tau_red=Tau(self.reduced.values)
        result = []
        if V.type == 'kron':
            for idx in itertools.product(*(range(di) for di in tau_red.G)): # Choice of one index in each component of tau ie a row in each column of the partial matrix
                if sum(tau_red.components[j][i] for j,i in enumerate(idx))  == 0:
                    result.append(idx)
            return result

        # In this case it a kind of mix between all-weights for boson and orthogonal weights
        s=len(self.reduced.values[0])
        for i,w in enumerate(itertools.combinations(range(s+V.nb_part-1),s-1)) :
                L=[w[0]]
                for j,p in enumerate(itertools.pairwise(w)):
                    L.append(p[1]-p[0]-1)
                #Lself.G.rank+V.nb_part-w[-1]-2
                L.append(s+V.nb_part-w[-1]-2)
                # check the weight condition
                if sum(i*t for i,t in zip(L,self.reduced.values[0])) == 0:
                    result.append(L)                        
        return result

class ReducedTau:
    """
    Tau in a reduction form
    
    Example:
    >>> tau = Tau(((3, 3, 2, 2), (2, 2, 1), (2, 2, 1),(1,)))
    >>> tau
    3 3 2 2 | 2 2 1 | 2 2 1 | 1
    >>> rtau = tau.reduced
    >>> rtau
    3^2 2^2 | 2^2 1^1 | 2^2 1^1 | 1^1
    >>> rtau.values
    3 2 | 2 1 | 2 1 | 1
    >>> rtau.mult
    2 2 | 2 1 | 2 1 | 1
    """
    __slots__ = 'values', 'mult'
    values: Blocks[int]
    mult: Blocks[int]

    def __init__(self, tau: Tau):
        from .utils import group_by_block
        values = Blocks.from_flatten([], (0,) * len(tau.G))
        mult = Blocks.from_flatten([], (0,) * len(tau.G))

        for j, component in enumerate(tau.components):
            values[j], mult[j] = zip(*group_by_block(component))

        self.values = values.freeze()
        self.mult = mult.freeze()

    def __len__(self) -> int:
        """ Number of components """
        return len(self.values)

    @property
    def small_d(self) -> LinGroup:
        return LinGroup(self.values.sizes)

    # FIXME: it seems that i,j is not really needed => simplification?
    def __getitem__(self, idx: tuple[int, int]) -> tuple[int, int]:
        return self.values[idx[1], idx[0]], self.mult[idx[1], idx[0]]

    def __repr__(self) -> str:
        return f" | ".join(
            " ".join(f"{v}^{m}" for v, m in zip(cv, cm))
            for cv, cm in zip(self.values.blocks, self.mult.blocks)
        )

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ReducedTau):
            return NotImplemented
        return self.values == other.values and self.mult == other.mult

    def __hash__(self) -> int:
        """ Hash consistent with equality so that to be safely used in a set or a dict """
        return hash((self.values, self.mult))

    
def unique_modulo_symmetry_list_of_tau(seq_tau: Iterable[Tau]) -> set[Tau]:
    """
    Unique sequence of tau modulo the it's symmetries

    Example:
    >>> G = LinGroup([2, 2, 2, 1, 1,1])
    >>> t1 = Tau.from_flatten([1, 6, 2, 1, 5, 1, 4, 5, 3], G)
    >>> t2 = t1.end0_representative
    >>> t3 = t2.sort_mod_sym_dim
    >>> t4 = t1.sort_mod_sym_dim
    >>> t5 = Tau.from_flatten([2, 6, 2, 1, 5, 1, 4, 5, 3], G)
    >>> for tau in unique_modulo_symmetry_list_of_tau((t1, t2, t3, t4, t5)):
    ...     print(tau)
    -4 0 | -3 0 | 4 0 | 0 | 0 | 21
    -4 0 | -3 0 | 4 0 | 0 | 0 | 20
    """
    return {tau.end0_representative.sort_mod_sym_dim for tau in seq_tau}

def full_under_symmetry_list_of_tau(seq_tau: Iterable[Tau]) -> Iterable[Tau]:
    """ provides all the elements in orbits under symmetry for elements in seq_tau
    """
    return itertools.chain.from_iterable(tau.orbit_symmetries() for tau in seq_tau)


def find_1PS(V: Representation, quiet: bool = False) -> Sequence["Tau"]:
    """
    Same as find_1PS_reg_mod_sym_dim without regularity condition
    Computed by 
    """
    
    # Initialisation with regular 1-PS
    List_1PS = []
    
    if V.type=='kron' :
        # List of representations corresponding to various tori S
        sub_rep=[Representation(LinGroup(list(p)),'kron') for p in Partition(list(V.G)).all_subpartitions()][1:] #[1:] excludes 1... 1
        for Vred in sub_rep:
            umax=V.G.u_max(Vred.G)
            #Recover by induction all candidates 1-PS mod symmetry
            List_1PS_Vred_reg=[]
            for H in find_hyperplanes_reg_mod_outer(Vred.all_weights, Vred, umax): # Important to keep sym at None
                taured=Tau.from_zero_weights(H, Vred)
                if taured.is_dom_reg : # We keep only dominant regular 1-PS
                    List_1PS_Vred_reg+=[t for t in taured.orbit_symmetries()]
                elif taured.opposite.is_dom_reg :
                    List_1PS_Vred_reg+=[t for t in taured.opposite.orbit_symmetries()]
                
            # Suppress the repetitions (one hyperplane can be spanned by several collections of weights)       
            List_1PS_Vred_reg=list(set(List_1PS_Vred_reg))
            
            #List_1PS_smalld_reg_mod_outer= list(find_1PS_reg_mod_sym_dim(Vred,umax))
            if not quiet:
                print('For G=',Vred.G,'we get',len(List_1PS_Vred_reg),' candidates regular dominant')

            #List_1PS_smalld_reg=sum([list(tau.orbit_symmetries()) for tau in List_1PS_smalld_reg_mod_sym]  ,[])
            List_1PS_Vred_extended=[]
            for permut in Permutation.embeddings_mod_sym(V.G, Vred.G):
                for tau in List_1PS_Vred_reg:
                    tau_twist=Tau([tau.components[i] for i in permut])
                    list_tau_extended=tau_twist.m_extend_with_repetitions(V.G)
                    for tau_ext in list_tau_extended:
                        if len(flatten_dictionary(tau_ext.positive_weights(V)))<=tau_ext.dim_Pu:
                            List_1PS_Vred_extended.append(tau_ext)
            List_1PS+=unique_modulo_symmetry_list_of_tau(List_1PS_Vred_extended)
            tau_ext=List_1PS[0]
    else :
        list_partS=[p for p in Partition.all_for_integer(V.G.rank)][1:] #[1:] excludes n, so S is the center of G
        for partS in list_partS :
            from math import floor
            umax=floor((V.G.dim-sum([x**2 for x in partS]))/2) # dim of P^u(tau)
            weights=[chi for chi in V.weights_of_S(partS)]
            Gred=LinGroup([len(partS)])
            Vred=Representation(Gred,V.type,V.nb_part)
            sym=list(symmetries(partS))
            for H in find_hyperplanes_reg_mod_outer(weights, Vred, umax, sym):
                taured=Tau.from_zero_weights(H, Vred)
                tau=taured.extend_from_S(partS)
                # dominant 1-PS corresponding to tau and -tau
                l1=list(tau.flattened)
                l1.sort(reverse=True)
                l2=list(tau.opposite.flattened)
                l2.sort(reverse=True)
                List_1PS+=[Tau.from_flatten(l1,V.G),Tau.from_flatten(l2,V.G)]

    List_1PS=list(set(List_1PS)) # Suppressing the repetitions
    
    # Checking that the candidates really give a candidate
        
    ListRes=[tau for tau in List_1PS if check_hyperplane_dim([chi for chi in tau.orthogonal_weights(V)], V.dim_cone-1)] # TODO : liste de chi utile ou itérable suffit ?
    
    return ListRes
