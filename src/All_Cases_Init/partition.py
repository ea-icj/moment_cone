from .typing import *
from .utils import is_decreasing, trim_zeros

import itertools

__all__ = (
    "OurPartition",
)

class OurPartition:
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
    def all_for_integer(n: int) -> Iterable["OurPartition"]: #TODO : est utilisé plusieurs fois pour le même n. staticmethod est adapté ?
        """
        Generates all partitions of an integer n >= 0.

        Could be optimized but it is clearly enough for the n we will consider.
        """
        if n <= 0:
            yield OurPartition((), check=False)
            return
    
        for head in range(n, 0, -1):
            for tail in OurPartition.all_for_integer(n - head):
                if len(tail) == 0 or head >= tail[0]:
                    yield OurPartition((head,) + tail._data, check=False)

    @staticmethod
    def all_of_height(height: int, lambda_max: int) -> Iterable["OurPartition"]:
        """
        Generates all partitions of given height and with given maximum value (included).
        """
        # Note that combinations_with_replacement keeps order of input sequence
        for w in itertools.combinations_with_replacement(reversed(range(lambda_max + 1)), height):
            yield OurPartition(w, check=False)

    # TODO: property? cached?
    def all_subpartitions(self) -> Iterable["OurPartition"]:
        """
        All sub-partitions of the partition
        
        Example:
        >>> p = OurPartition((4, 2, 2, 1))
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
            yield OurPartition((), check=False)

        for x in range(self[0]):
            tail = OurPartition([min(x + 1, y) for y in self._data[1:]])
            for tail_sp in tail.all_subpartitions():
                yield OurPartition((x + 1, *tail_sp), check=False)

    
