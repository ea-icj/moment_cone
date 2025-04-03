__all__ = (
    "Root",
)

from dataclasses import dataclass
import itertools

from .rings import Vector
from .linear_group import LinearGroup
from .typing import *


@dataclass(frozen=True, slots=True)
class Root:
    """ Root element for tau """
    k: int
    i: int
    j: int

    def __hash__(self) -> int:
        return hash((self.k, self.i, self.j))
    
    def __lt__(self, other):
        """Lexicographic order on (i, j, k)."""
        return (self.i, self.j, self.k) < (other.i, other.j, other.k)


    @property
    def is_in_U(self) -> bool:
        """ Check if this root is in U """
        return self.i < self.j
    
    @property
    def opposite(self) -> "Root":
        """ Return the opposite of a root (k,i,j -> kj,i) """
        return Root(self.k, self.j, self.i)

    def to_vector(self, G: LinearGroup) -> Vector:
        """
        Returns self as a vector in Z**(G.rank). A kind of flatten.
        
        """
        from .rings import vector, QQ
        v = vector(QQ, G.rank)
        shift=sum(G[:self.k])
        v[shift+self.i]=1
        v[shift+self.j]=-1
        return(v)
    
    @staticmethod
    def all_of_U(G: LinearGroup) -> Iterable["Root"]:
        """
        Returns all possible root from U for given dimensions
        
        Example:
        >>> G = LinearGroup((2, 3))
        >>> for root in Root.all_of_U(G):
        ...     print(root)
        Root(k=0, i=0, j=1)
        Root(k=1, i=0, j=1)
        Root(k=1, i=0, j=2)
        Root(k=1, i=1, j=2)
        """
        for k, dk in enumerate(G):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)

    @staticmethod
    def all_of_B(G: LinearGroup) -> Iterable["Root"]:
        """
        Returns all possible root from B for a given group
        
        FIXME: verify example

        Example:
        >>> G = LinearGroup((2, 3))
        >>> for root in Root.all_of_B(G):
        ...     print(root)
        Root(k=0, i=0, j=0)
        Root(k=0, i=0, j=1)
        Root(k=0, i=1, j=1)
        Root(k=1, i=0, j=0)
        Root(k=1, i=0, j=1)
        Root(k=1, i=0, j=2)
        Root(k=1, i=1, j=1)
        Root(k=1, i=1, j=2)
        Root(k=1, i=2, j=2)
        """
        for k, dk in enumerate(G):
            for i in range(dk):
                for j in range(i, dk):
                    yield Root(k, i, j)
                
    @staticmethod
    def all_of_K(G: LinearGroup) -> Iterable["Root"]:
        """
        Returns all possible roots for Lie(K) for given dimensions. Root(k,i,i) allowed corresponding to diagonal matrices. 
        We get a set indexing a bases of Lie(K)
        """
        for k,e in enumerate(G):
            for i,j in itertools.product(range(e), repeat=2):
                yield Root(k,i,j) 

    def index_in_all_of_K(self, G: LinearGroup, use_internal_index: bool = True) -> int:
        """
        Returns index of this weight in the lexicographical order for given dimensions (see `all_of_K` method)
        
        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """

        tot=sum(x**2 for x in G[:self.k])
        return tot+G[self.k]*self.i+self.j
    
    @staticmethod
    def all(G: LinearGroup) -> Iterable["Root"]:
        """
        Returns all possible roots from G (i != j) for given dimensions
        
        Example:
        >>> G = LinearGroup((2, 3))
        >>> for root in Root.all(G):
        ...     print(root)
        Root(k=0, i=0, j=1)
        Root(k=0, i=1, j=0)
        Root(k=1, i=0, j=1)
        Root(k=1, i=1, j=0)
        Root(k=1, i=0, j=2)
        Root(k=1, i=2, j=0)
        Root(k=1, i=1, j=2)
        Root(k=1, i=2, j=1)
        """
        for k, dk in enumerate(G):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)
                yield Root(k, j, i)


    @staticmethod
    def all_of_T(G: LinearGroup) -> Iterable["Root"]:
        """
        Returns all possible roots from T (i == j) for given dimensions
        
        Example:
        >>> G = LinearGroup((2, 3))
        >>> for root in Root.all_of_T(G):
        ...     print(root)
        Root(k=0, i=0, j=0)
        Root(k=0, i=1, j=1)
        Root(k=1, i=0, j=0)
        Root(k=1, i=1, j=1)
        Root(k=1, i=2, j=2)
        """
        for k, dk in enumerate(G):
            for i in range(dk):
                yield Root(k, i, i)       
    
    @property
    #TODO: erase this one if not used in Groebner
    def short_rep(self) -> tuple[int, int, int]:
        return self.k, self.i, self.j
    
    @property
    def as_list(self) -> list[int]:
        return [self.k,self.i,self.j]
