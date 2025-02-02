from .typing import *
from .utils import is_decreasing, trim_zeros

import itertools

__all__ = (
    "Partition", "ListPartPlus",
)

class Partition:
    """ Decreasing sequence of positive int """
    __slots__ = "_data",
    _data: tuple[int, ...]

    def __init__(self, p: Sequence[int], check: bool = True):
        # Auto trim the partition
        self._data = tuple(trim_zeros(p))
        assert not check or self.is_valid, "Invalid partition"

    @property
    def is_valid(self) -> bool:
        """ Check if partition is decreasing and positive """
        return len(self) == 0 or (is_decreasing(self._data) and self._data[-1] >= 0)

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, idx: int):
            return self._data[idx] if idx < len(self) else 0
        
    def __iter__(self) -> Iterator[int]:
        return iter(self._data)
    
    def __repr__(self) -> str:
        return f"Partition({self._data})"
    
    def __eq__(self, other) -> bool:
        return len(self) == len(other) and self._data == other._data
    
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

    def lambda_check(self,l : int) -> "Partition" :
        """
        l for GL(l). Max length of la. TODO : ajouter un assert
        """
        x = self[0]
        return(Partition([x-self[i] for i in range(l-1, 0, -1)]))

    def lambda_red(self,l : int) -> "Partition" : 
        """
        l for GL(l). Tensor with det to reduce la. The ouput can be thought as a representation of SL(l)
        """
        x = self[l-1]
        return(Partition([self[i]-x for i in range(l-1)]))


        

    
class ListPartPlus: 
    """
    A list of partitions with two more properties 
    - indices (a list of integers) and 
    - mult (an integer which is a multiplicity in representation theory.
    """

    def __init__(self, L : list[Partition], c : int,indices : list[int]=None):
        """
        Initializes an instance of ListPartPlus.
        """
        self.parts=L
        self.mult=c
        self.indices=indices
        

    def __repr__(self) -> str:
        if self.indices != None :
            return 'Partitions: '+str(self.parts)+', Indices: '+str(self.indices)+', Multiplicity: '+str(self.mult)
        else :
            return 'Partitions: '+str(self.parts)+', Multiplicity: '+ str(self.mult)

    def __eq__(self,other : "ListPartPlus") -> bool:
        if all(l==m for l,m in zip(self.parts,other.parts)) and self.mult==other.mult and self.indices==other.indices :
            return True
        return False
