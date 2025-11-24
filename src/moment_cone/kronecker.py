"""
Tools to compute the Kronecker coefficient of multiple partitions in an
efficient way using a cache
"""

__all__ = (
    "kronecker_product",
    "KroneckerCoefficient",
    "KroneckerCoefficientCache",
    "KroneckerCoefficientMLCache",
)

from sage.combinat.sf.sf import SymmetricFunctions # type: ignore
from sage.rings.rational_field import QQ # type: ignore

from .typing import *
from .partition import Partition


# Symmetric Functions over Rational Field in the Schur basis
__sym_f = SymmetricFunctions(QQ).schur()

def sym_f(p: Partition) -> Any:
    return __sym_f(tuple(p))


def kronecker_product(partitions: Sequence[Partition]) -> dict[Partition, int]:
    """
    Returns the Kronecker product decomposition of the given partitions

    The output is a dictionary of coefficient for each partition
    of the decomposition.

    >>> a = Partition(3, 2, 1)
    >>> b = Partition(2, 2, 2)
    >>> c = Partition(2, 2, 1, 1)
    >>> product = kronecker_product((a, b, c))
    >>> for p, coeff in sorted(product.items()):
    ...     print(f"{p}: {coeff}")
    Partition((6,)): 1
    Partition((5, 1)): 5
    Partition((4, 2)): 9
    Partition((4, 1, 1)): 10
    Partition((3, 3)): 5
    Partition((3, 2, 1)): 16
    Partition((3, 1, 1, 1)): 10
    Partition((2, 2, 2)): 5
    Partition((2, 2, 1, 1)): 9
    Partition((2, 1, 1, 1, 1)): 5
    Partition((1, 1, 1, 1, 1, 1)): 1

    >>> kronecker_product((Partition(3, 2, 1),))
    {Partition((3, 2, 1)): 1}

    >>> kronecker_product(())
    {}
    """
    if len(partitions) == 0:
        return dict()
    
    # Kronecker product
    product = sym_f(partitions[0])
    for p in partitions[1:]:
        product = product.kronecker_product(sym_f(p))

    # Splitting decomposition
    result: dict[Partition, int] = dict()
    for monomial, coeff in product.monomial_coefficients().items():
        # Forcing int coefficients to remove Sage custom Integer type
        result[Partition(map(int, monomial))] = int(coeff)

    return result


class KroneckerCoefficient:
    """
    Base class for computing the Kronecker coefficient of partitions.

    >>> kc = KroneckerCoefficient()
    >>> partitions = (3, 1, 1, 1, 1, 1), (4, 3, 1), (2, 2, 2, 2), (5, 1, 1, 1), (4, 2, 1,1)
    >>> partitions = tuple(Partition(p) for p in partitions)
    >>> kc(partitions)
    1607
    >>> kc(partitions[:3])
    1
    >>> kc(partitions[:2])
    0
    >>> kc(partitions[:1])
    1
    >>> kc((partitions[0], partitions[0]))
    1
    """
    def __call__(self, partitions: Sequence[Partition]) -> int:
        """ Returns the Kronecker coefficient of given partitions
        
        May rely on symmetries and other properties.
        """
        if len(partitions) == 0:
            return 0
        elif len(partitions) == 1:
            return 1
        elif len(partitions) == 2:
            return partitions[0] == partitions[1]
        else:
            return self._kernel(self._sort(partitions))

    def product(self, partitions: Sequence[Partition]) -> dict[Partition, int]:
        """ Computes the Kronecker product of given partitions """
        if len(partitions) == 0:
            return dict()
        elif len(partitions) == 1:
            return {partitions[0]: 1}
        else:
            return self._product(self._sort(partitions))
    
    def _sort(self, partitions: Iterable[Partition]) -> list[Partition]:
        """
        Internal method to sort partitions so that to reduce requests
        
        Ideally, to reduce the size of the cache, a same set of partitions but
        in different order should lead to the same ordered list returned by
        this method.

        Sorting by increasing length and lexicographical order is faster
        and requires a smaller cache.
        """
        return sorted(partitions, key=lambda p: (len(p), p))
    
    def _product(self, partitions: Sequence[Partition]) -> dict[Partition, int]:
        """ Overridable kernel that computes the Kronecker product of given partitions """
        return kronecker_product(partitions)
    
    def _kernel(self, partitions: Sequence[Partition]) -> int:
        """ Overridable kernel that computes the Kronecker coefficient """
        return self._product(partitions[:-1]).get(partitions[-1], 0)

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


class KroneckerCoefficientCache(KroneckerCoefficient):
    """
    Kronecker coefficient computation using cache to store the computed decompositions

    >>> kc = KroneckerCoefficientCache()
    >>> print(kc)
    KroneckerCoefficientCache(#cache=0 (), #hit=0, #miss=0)

    >>> partitions = (3, 1, 1, 1, 1, 1), (4, 3, 1), (2, 2, 2, 2), (5, 1, 1, 1), (4, 2, 1, 1)
    >>> partitions = tuple(Partition(p) for p in partitions)
    >>> kc(partitions)
    1607
    >>> print(kc)
    KroneckerCoefficientCache(#cache=1 (#4=1), #hit=0, #miss=1)

    >>> kc((partitions[0], partitions[1], partitions[3]))
    1
    >>> print(kc)
    KroneckerCoefficientCache(#cache=2 (#2=1,#4=1), #hit=0, #miss=2)

    >>> kc(partitions[:3])
    1
    >>> kc(partitions[:2])
    0
    >>> kc(partitions[:1])
    1
    >>> kc((partitions[0], partitions[0]))
    1
    """
    _cache: dict[tuple[Partition, ...], dict[Partition, int]]
    _hit: int
    _miss: int

    def __init__(self, file_prefix: Optional[str] = None):
        super().__init__()
        self._cache = dict()
        self._hit = 0
        self._miss = 0

        if file_prefix is not None:
            self.load_cache(file_prefix)

    def _product(self, partitions: Sequence[Partition]) -> dict[Partition, int]:
        partitions = tuple(partitions)
        try:
            product = self._cache[partitions]
            self._hit += 1
            return product
        except KeyError:
            self._miss += 1
            product = kronecker_product(partitions)
            self._cache[partitions] = product
            return product

    def __repr__(self) -> str:
        from .utils import group_by_block
        cache_details = group_by_block(sorted(len(key) for key in self._cache.keys()))
        cache_details_str = ",".join(f"#{n}={c}" for n, c in cache_details)
        return (
            f"{type(self).__name__}("
            f"#cache={len(self._cache)} ({cache_details_str}), #hit={self._hit}, #miss={self._miss}"
            ")"
        )
    
    def load_cache(self, file_prefix: str, clear: bool = False) -> None:
        """ Load cache from given filename prefix (will append .pkl.xz) """
        if clear:
            self._cache = dict()

        import lzma, pickle, re
        file_name = file_prefix + ".pkl.xz"
        if re.match(r"^https?://", file_name):
            import urllib.request
            with urllib.request.urlopen(file_name) as rfh:
                with lzma.open(rfh, "rb") as fh:
                    self._cache.update(pickle.load(fh))
        else:
            with lzma.open(file_name, "rb") as fh:
                self._cache.update(pickle.load(fh))

    def save_cache(self, file_prefix: str) -> None:
        """ Save cache in given filename prefix (will append .pkl.xz) """
        import lzma, pickle
        with lzma.open(file_prefix + ".pkl.xz", "wb") as fh:
            pickle.dump(self._cache, fh)


class KroneckerCoefficientMLCache(KroneckerCoefficientCache):
    """
    Kronecker coefficient computation using a multi-level cache.
    
    The cache contains the decompositions of the full Kronecker product of
    encountered sequence of partitions and also the decompositions of all
    intermediary product.

    It uses a divide-to-conquer strategy to reduce the length of the products
    that needed to be cached.

    Furthermore, the last step that compute the Kronecker coefficient is
    optimized from the decompositions of the sequence of partitions splitted
    in two parts.

    >>> kc = KroneckerCoefficientMLCache()
    >>> print(kc)
    KroneckerCoefficientMLCache(#cache=0 (), #hit=0, #miss=0)

    >>> partitions = (3, 1, 1, 1, 1, 1), (4, 3, 1), (2, 2, 2, 2), (5, 1, 1, 1), (4, 2, 1,1)
    >>> partitions = tuple(Partition(p) for p in partitions)
    >>> kc(partitions)
    1607
    >>> print(kc)
    KroneckerCoefficientMLCache(#cache=3 (#2=2,#3=1), #hit=0, #miss=3)

    >>> kc((partitions[0], partitions[1], partitions[3]))
    1
    >>> print(kc)
    KroneckerCoefficientMLCache(#cache=3 (#2=2,#3=1), #hit=1, #miss=3)

    >>> kc(partitions[:3])
    1
    >>> kc(partitions[:2])
    0
    >>> kc(partitions[:1])
    1
    >>> kc((partitions[0], partitions[0]))
    1
    """
    def _product(self, partitions: Sequence[Partition]) -> dict[Partition, int]:
        partitions = tuple(partitions)
        result: dict[Partition, int]

        try:
            result = self._cache[partitions]
            self._hit += 1
            return result
        except KeyError:
            self._miss += 1
            head: Any
            tail: Any
            if len(partitions) == 2:
                head = sym_f(partitions[0])
                tail = sym_f(partitions[1])
            else:
                # Split sequence of partitions in two parts of close length
                # with the constraint that the first part must contains at least
                # two partitions.

                # Decomposition of the first part as a Sage object
                head_length = max(2, len(partitions) // 2)
                head = sum(
                    coeff * sym_f(p)
                    for p, coeff in self._product(partitions[:head_length]).items()
                )

                # Decomposition of the second part as a Sage object
                tail_length = len(partitions) - head_length
                if tail_length == 1:
                    tail = sym_f(partitions[-1])
                else:
                    tail = sum(
                        coeff * sym_f(p)
                        for p, coeff in self._product(partitions[head_length:]).items()
                     )

            product = head.kronecker_product(tail)

            # Converting decomposition to the Partition -> coeff dictionary
            result = dict()
            for monomial, coeff in product.monomial_coefficients().items():
                result[Partition(tuple(map(int, monomial)))] = int(coeff)

            self._cache[tuple(partitions)] = result
            return result
        
    def _kernel(self, partitions: Sequence[Partition]) -> int:
        if len(partitions) == 2:
            return 1 if partitions[0] == partitions[1] else 0
        
        if len(partitions) == 3:
            return self._product(partitions[:2]).get(partitions[2], 0)
        
        # For length at most 4, we can split the sequence in two parts
        # of at least two partitions each, compute the decomposition of each
        # part and then compute the Kronecker coefficient using a kind
        # of scalar product of the two decompositions.
        head_length = max(2, len(partitions) // 2)
        head = self._product(partitions[:head_length])
        tail = self._product(partitions[head_length:])
        return sum(head[m] * tail[m] for m in head.keys() & tail.keys())    
