__all__ = (
    "Partition",
    "EnhancedPartitionList",
    "gen_partitions"
)

import itertools

from .typing import *
from .utils import is_decreasing, trim_zeros

class Partition:
    """
    Decreasing sequence of positive int
    
    Example:
    >>> p = Partition((3, 2, 1))
    >>> print(p)
    Partition((3, 2, 1))

    >>> Partition(2, 2, 2) # Also possible without inner parenthesis
    Partition((2, 2, 2))
    """
    __slots__ = "_data",
    _data: tuple[int, ...]

    def __init__(self, p: int | Iterable[int], *tail: int, check: bool = True):
        if isinstance(p, Iterable):
            assert len(tail) == 0
            coeffs = tuple(p)
        else:
            coeffs = (p,) + tail
        
        # Auto trim the partition
        self._data = cast(tuple[int, ...], trim_zeros(coeffs))
        assert not check or self.is_valid, "Invalid partition"

    @property
    def is_valid(self) -> bool:
        """ Check if partition is decreasing and positive """
        return len(self) == 0 or (is_decreasing(self._data) and self._data[-1] >= 0)

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, idx: int) -> int:
            return self._data[idx] if idx < len(self) else 0
        
    def __iter__(self) -> Iterator[int]:
        return iter(self._data)
    
    def __repr__(self) -> str:
        return f"Partition({self._data})"
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Partition):
            return NotImplemented
        return self._data == other._data
    
    def __lt__(self, other: object) -> bool:
        """ Reverse lexicographical order """
        if not isinstance(other, Partition):
            return NotImplemented
        return self._data > other._data

    def __le__(self, other: object) -> bool:
        """ Reverse lexicographical order """
        if not isinstance(other, Partition):
            return NotImplemented
        return self._data >= other._data
    
    def __add__(self, other: object) -> "Partition":
        """ Concatenating a Partition with other iterable
        
        Example:
        >>> p = Partition((3, 2))
        >>> p + Partition((2, 1, 1))
        Partition((3, 2, 2, 1, 1))
        >>> p + [2, 1, 1]
        Partition((3, 2, 2, 1, 1))
        """
        if not isinstance(other, Iterable):
            return NotImplemented
        return Partition(itertools.chain(self, other))

    def __hash__(self) -> int:
        return hash(self._data)

    def pad(self, length: int) -> tuple[int, ...]:
        """ Returns a padded version of this partition """
        assert length >= len(self._data), "Padding length must be greater that Partition length"
        return self._data + (0,) * (length - len(self._data))

    @staticmethod
    def all_for_integer(n: int, max_length: Optional[int] = None) -> Iterable["Partition"]:
        """
        Generates all partitions of an integer n >= 0 with optional maximal length.

        Could be optimized but it is clearly enough for the n we will consider.
        """        
        if n <= 0:
            yield Partition((), check=False)
            return

        if max_length is None:
            max_length = n
        if max_length == 0:
            # If max_length is zero and n is still strictly greater than 0,
            # then it means that the current branch don't lead to valid Partitions
            return

        for head in range(n, 0, -1):
            for tail in Partition.all_for_integer(n - head, max_length - 1):
                if len(tail) == 0 or head >= tail[0]:
                    yield Partition((head,) + tail._data, check=False)

    @staticmethod
    def all_of_length(length: int, lambda_max: int) -> Iterable["Partition"]:
        """
        Generates all partitions of given length and with given maximum value (included).
        """
        # Note that combinations_with_replacement keeps order of input sequence
        for w in itertools.combinations_with_replacement(reversed(range(lambda_max + 1)), length):
            yield Partition(w, check=False)

    # TODO: property? cached?
    def all_subpartitions(self) -> Iterable["Partition"]:
        """
        All sub-partitions of the partition
        
        Example:
        >>> p = Partition((4, 2, 2, 1))
        >>> for sp in p.all_subpartitions():
        ...     print(sp)
        Partition((1, 1, 1, 1))
        Partition((2, 1, 1, 1))
        Partition((2, 2, 1, 1))
        Partition((2, 2, 2, 1))
        Partition((3, 1, 1, 1))
        Partition((3, 2, 1, 1))
        Partition((3, 2, 2, 1))
        Partition((4, 1, 1, 1))
        Partition((4, 2, 1, 1))
        Partition((4, 2, 2, 1))
        """
        if len(self) == 0 :
            yield Partition((), check=False)

        for x in range(self[0]):
            tail = Partition([min(x + 1, y) for y in self._data[1:]])
            for tail_sp in tail.all_subpartitions():
                yield Partition((x + 1, *tail_sp), check=False)
        

    def lambda_check(self, l: int) -> "Partition":
        """
        l for GL(l). Max length of la. TODO : ajouter un assert
        """
        x = self[0]
        return Partition([x - self[i] for i in range(l - 1, 0, -1)])

    def lambda_red(self, l: int) -> "Partition": 
        """
        l for GL(l). Tensor with det to reduce la. The output can be thought as a representation of SL(l)
        """
        if l == 0:
            return(self)
        x = self[l - 1]
        return Partition([self[i] - x for i in range(l - 1)])

    @staticmethod
    def join(partitions: "Iterable[Partition]") -> "Partition":
        """
        Returns the component-wise maximum of the given partitions
        
        Examples:
        >>> Partition.join((
        ...     Partition(3, 2, 2),
        ...     Partition(4, 1, 1),
        ...     Partition(3, 3, 2),
        ... ))
        Partition((4, 3, 2))
        
        >>> Partition.join((
        ...     Partition(3, 2, 2, 2, 1),
        ...     Partition(4, 1, 1, 1),
        ...     Partition(3, 3, 3, 1, 1, 1),
        ... ))
        Partition((4, 3, 3, 2, 1, 1))
        """
        return Partition(
            map(max, itertools.zip_longest(*partitions, fillvalue=0))
        )


def find_max_index(L: Partition, Out: Partition) -> Optional[int]:
    """
    Biggest index i that can be increased in L to stay a sub-partition of Out by keeping the first values of L unchanged.
    Return None if such an index does not exist.

    Example:

    >>> find_max_index([2,1,1], [3,2,1])
    1
    """
    for i in reversed(range(1, len(Out))):  # Parcours à l'envers
        if L[i - 1] > L[i] and Out[i] > L[i]:
            return i
    if Out[0] > L[0] :
        return 0
    return None


def gen_partitions(
        min_weight: int,
        max_weight: int,
        In: Partition,
        Out: Partition
    ) -> list[Partition]:
    """
    Generate all the partitions containing In contained in Out and of weight in [min_weight, max_weight]
    >>> for l in gen_partitions(4, 6, [2,1,0], [4,3,1]):
    ...     print(l)
    Partition((2, 1, 1))
    Partition((2, 2))
    Partition((2, 2, 1))
    Partition((3, 1))
    Partition((3, 1, 1))
    Partition((3, 2))
    Partition((3, 2, 1))
    Partition((3, 3))
    Partition((4, 1))
    Partition((4, 1, 1))
    Partition((4, 2))
    """
    if len(Out)==0:
        if min_weight == 0 :
            return [Partition([0])]
        else :
            return []
    result = []    
    current = [x for x in In]+[0]*(len(Out)-len(In)) # copy
    S=sum(current)
    if S<=max_weight and S>=min_weight :
        result.append(Partition([x for x in current]))
    i = find_max_index(Partition(current), Out)
    while i is not None:
        current[i]+=1
        for j in range(i+1,len(Out)):
            current[j]=In[j]
        i=find_max_index(Partition(current), Out)
        S=sum(current)
        if S<=max_weight and S>=min_weight :
            result.append(Partition([x for x in current]))
    return result    

class EnhancedPartitionList:
    """
    A list of partitions additional properties.
    
    - indices (a list of integers) and 
    - mult (an integer which is a multiplicity in representation theory.
    """
    partitions: list[Partition]
    mult: int
    indices: Optional[list[int]]

    def __init__(
            self,
            partitions: Iterable[Partition],
            mult: int,
            indices: Optional[Iterable[int]] = None):
        """
        Initializes an instance of ListPartPlus.
        """
        self.partitions = list(partitions)
        self.mult = mult
        if indices is None:
            self.indices = None
        else:
            self.indices = list(indices)        

    def __repr__(self) -> str:
        if self.indices != None :
            return 'Partitions: '+str(self.partitions)+', Indices: '+str(self.indices)+', Multiplicity: '+str(self.mult)
        else :
            return 'Partitions: '+str(self.partitions)+', Multiplicity: '+ str(self.mult)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, EnhancedPartitionList):
            return NotImplemented
        
        return self.partitions == other.partitions and self.mult==other.mult and self.indices==other.indices
