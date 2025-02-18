from .typing import *
from .dimension import Dimension
from .rings import QQ, vector, Vector

from dataclasses import dataclass
import itertools

__all__ = (
    "Root",
)

@dataclass(frozen=True, slots=True)
class Root:
    """ Root element for tau """
    k: int
    i: int
    j: int

    def __hash__(self) -> int:
        return hash((self.k, self.i, self.j))
    
    @property
    def is_in_U(self) -> bool:
        """ Check if this root is in U """
        return self.i < self.j
    
    
    @property
    def short_repr(self) -> Iterable[int]:
        """  returns a short representation for our root """
        return [self.k,self.i,self.j]
    
    @property
    def opposite(self) -> "Root":
        """ Return the opposite of a root (k,i,j -> kj,i) """
        return Root(self.k, self.j, self.i)

    
    def to_vector(self, d: Dimension) -> Vector:
        """
        Returns self as a vector in Z**(sum(d) + 1.
        
        A kind of flatten.

        Example:
        >>> d = Dimension((2, 3))
        >>> root = Root(1, 0, 2)
        >>> root.to_vector(d)
        (0, 0, 0, 1, 0, -1)
        """
        v = vector(QQ, d.sum + 1)
        shift = 1+sum(d[:self.k]) # 1 because the fist term is for C^* and has not root
        v[shift + self.i] = 1
        v[shift + self.j] = -1
        return v
    
    @staticmethod
    def all_of_U(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible root from U for given dimensions
        
        Example:
        >>> d = Dimension((2, 3))
        >>> for root in Root.all_of_U(d):
        ...     print(root)
        Root(k=0, i=0, j=1)
        Root(k=1, i=0, j=1)
        Root(k=1, i=0, j=2)
        Root(k=1, i=1, j=2)
        """
        for k, dk in enumerate(d):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)
                
    @staticmethod
    def all_of_K(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible roots for Lie(K) for given dimensions. Root(k,i,i) allowed corresponding to diagonal matrices. 
        We get a set indexing a bases of Lie(K)
        """
        for k,e in enumerate(d):
            for i,j in itertools.product(range(e), repeat=2):
                yield Root(k,i,j) 

    def index_in_all_of_K(self, d: Dimension, use_internal_index: bool = True) -> int:
        """
        Returns index of this weight in the lexicographical order for given dimensions (see `all_of_K` method)
        
        By default, it will returns the index attribute (if not None) assuming that it has been defined
        for the same dimensions. `use_internal_index` can be set to `False` in order to force the computation
        of the index for the given dimension. In that case, the internal index will be updated for later reuse.
        """

        tot=sum(x**2 for x in d[:self.k])
        return tot+d[self.k]*self.i+self.j
    
    @staticmethod
    def all(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible roots from G (i != j) for given dimensions
        
        Example:
        >>> d = Dimension((2, 3))
        >>> for root in Root.all(d):
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
        for k, dk in enumerate(d):
            for i, j in itertools.combinations(range(dk), 2):
                yield Root(k, i, j)
                yield Root(k, j, i)


    @staticmethod
    def all_of_T(d: Dimension) -> Iterable["Root"]:
        """
        Returns all possible roots from T (i == j) for given dimensions
        
        
        Example:
        >>> d = Dimension((2, 3))
        >>> for root in Root.all_of_T(d):
        ...     print(root)
        Root(k=0, i=0, j=0)
        Root(k=0, i=1, j=1)
        Root(k=1, i=0, j=0)
        Root(k=1, i=1, j=1)
        Root(k=1, i=2, j=2)
        """
        for k, dk in enumerate(d):
            for i in range(dk):
                yield Root(k, i, i)            

