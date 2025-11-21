import itertools

from .typing import *

__all__ = (
    "Blocks",
)

class Blocks(Generic[T]):
    """
    Implementation of blocks of data of different sizes.

    Can be initialized for the list of the blocks (`Blocks.from_blocks`)
    or from a sequence of values and the sequence of size of each blocks
    (`Blocks.from_flatten` or using the default constructor).

    It is internally stored as a flattened sequence of values and the indexes of each blocks.

    When initializing from a sequence, the mutability of the sequence (ie list of tuple) defined the mutability of the instance of Blocks. When mutable, a whole block or an element of block can be modified, eg:

    >>> b = Blocks([1, 2, 3], (2, 1))
    >>> print(b)
    1 2 | 3
    >>> b[1] = [3, 4, 5]
    >>> print(b)
    1 2 | 3 4 5
    >>> b[0, 1] = 3
    >>> print(b)
    1 3 | 3 4 5

    When immutable, hash is consistent with equality thus allowing to safely use a Blocks in a set of a dict. An mutable Blocks can be converted to a immutable one using the `freeze` method:
    
    >>> b1 = Blocks([1, 2, 3], (2, 1))
    >>> b2 = Blocks([1, 2, 3], (2, 1))
    >>> b1 == b2, hash(b1), hash(b2) # doctest: +SKIP
    (True, 8789978014300, 8789977956485)
    >>> b1f = b1.freeze()
    >>> b2f = b2.freeze()
    >>> b1f == b2f, hash(b1f), hash(b2f)
    (True, 1557127941034686738, 1557127941034686738)
    """
    __slots__ = 'flatten', '_indexes'
    flatten: Sequence[T] #: Flattened version of the data (can be mutable or not)
    _indexes: tuple[int, ...] #: Boundaries of each block

    def __init__(self, flatten: Sequence[T], sizes: Iterable[int]):
        """
        Initialization from a flattened sequence and block sizes

        The resulting Blocks instance can be mutable if flatten is mutable.
        """
        self.flatten = flatten
        self._indexes = tuple(itertools.accumulate(sizes, initial=0))
        assert self._indexes[-1] == len(self.flatten), "Invalid sizes"

    @staticmethod
    def from_flatten(flatten: Sequence[T], sizes: Iterable[int]) -> "Blocks[T]":
        """
        Initialization from a flattened sequence and block sizes (same as __init__)
        
        The resulting Blocks instance can be mutable if flatten is mutable.
        """
        return Blocks(flatten, sizes)
    
    @staticmethod
    def from_blocks(blocks: Iterable[Sequence[T]]) -> "Blocks[T]":
        """ Initialization from multiple blocks (immutable) """
        blocks1, blocks2 = itertools.tee(blocks, 2) # To allow multiple reads of a generator
        return Blocks.from_flatten(
            tuple(itertools.chain.from_iterable(blocks1)),
            map(len, blocks2)
        )

    @property
    def is_frozen(self) -> bool:
        """ Return True if this instance is immutable """
        return not isinstance(self.flatten, MutableSequence)
    
    def freeze(self) -> "Blocks[T]":
        """ Return a freezed (immutable) version of this Blocks """
        if self.is_frozen:
            return self
        else:
            return Blocks.from_flatten(tuple(self.flatten), self.sizes)

    @property
    def sizes(self) -> tuple[int, ...]:
        """ Size of each block """
        return tuple(b - a for a, b in itertools.pairwise(self._indexes))
    
    @property
    def blocks(self) -> Iterable[Sequence[T]]:
        """ Tuple of blocks """
        return (self.flatten[a:b] for a, b in itertools.pairwise(self._indexes))

    def __len__(self) -> int:
        """ Number of blocks """
        return len(self._indexes) - 1
    
    def __iter__(self) -> Iterator[Sequence[T]]:
        """ Iterating through the blocks """
        return iter(self.blocks)
    
    def __repr__(self) -> str:
        """ String representation """
        return " | ".join(" ".join(map(str, c)) for c in self.blocks)
        
    @overload
    def __getitem__(self, idx: int) -> Sequence[T]:
        """ Get a block """
        ...

    @overload
    def __getitem__(self, idx: slice) -> Iterable[Sequence[T]]:
        """ Get a slice of blocks """
        ...

    @overload
    def __getitem__(self, idx: tuple[int, int]) -> T:
        """ Get an element of a block by its index (block_id, element_id) """
        ...

    def __getitem__(self, idx: int | slice | tuple[int, int]) -> Sequence[T] | Iterable[Sequence[T]] | T:
        """ Get a block or an element """
        # Get a block by index
        if not isinstance(idx, (slice, tuple)): # To take into account integral types that are not int (like Integer from Sage)
            if not 0 <= idx < len(self):
                raise IndexError("block index out of range")
            a, b = self._indexes[idx], self._indexes[idx + 1]
            return self.flatten[a:b]
        
        # Get a slice of blocks
        elif isinstance(idx, slice):
            start = idx.start or 0
            stop = idx.stop or len(self)
            step = idx.step or 1
            return (self[i] for i in range(start, stop, step))

        # Get an element of a block by its index (block_id, element_id)
        else:
            block_id, element_id = idx
            if not 0 <= block_id < len(self):
                raise IndexError("block index out of range")

            global_id = self._indexes[block_id] + element_id
            if element_id < 0 or global_id >= self._indexes[block_id + 1]:
                raise IndexError("element index out of range")

            return self.flatten[global_id]    
    
    @overload
    def __setitem__(self, idx: int, value: Sequence[T]) -> None:
        """ Set a block """
        ...

    @overload
    def __setitem__(self, idx: tuple[int, int], value: T) -> None:
        """ Set an element of a block by its index (block_id, element_id) """
        ...

    def __setitem__(self, idx: int | tuple[int, int], value: Sequence[T] | T) -> None:
        """ Set a block or an element """
        assert isinstance(self.flatten, MutableSequence), "this Blocks is frozen"

        # Modify a block
        if not isinstance(idx, (slice, tuple)): # To take into account integral types that are not int (like Integer from Sage)
            assert isinstance(value, Sequence)
            if not 0 <= idx < len(self):
                raise IndexError("block index out of range")
            a, b = self._indexes[idx], self._indexes[idx + 1]
            self.flatten[a:b] = value

            # Updating indexes
            shift = len(value) - (b - a)
            self._indexes = self._indexes[:idx + 1] + tuple(i + shift for i in self._indexes[idx + 1:])

        # Modify an element of a block by its index (block_id, element_id)        
        else:
            assert not isinstance(value, Sequence)
            block_id, element_id = idx
            if not 0 <= block_id < len(self):
                raise IndexError("block index out of range")

            global_id = self._indexes[block_id] + element_id
            if element_id < 0 or global_id >= self._indexes[block_id + 1]:
                raise IndexError("element index out of range")

            self.flatten[global_id] = value

    def __eq__(self, other: Any) -> bool:
        """
        Equality test
        
        Currently returns False if internal storage are not comparable (list and tuple)
        """
        if not isinstance(other, Blocks):
            return NotImplemented
        return self.flatten == other.flatten and self._indexes == other._indexes
    
    def __hash__(self) -> int:
        """ Hashing so that to be consistant with equality when immutable """
        if self.is_frozen:
            return hash((self.flatten, self._indexes))
        else:
            return super().__hash__()
    
