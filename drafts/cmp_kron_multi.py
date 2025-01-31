#!/usr/bin/env python3
from typing import Sequence, TypeAlias, Iterable, Iterator, Optional
from abc import ABC, abstractmethod
import itertools

from sage.all import SymmetricFunctions, QQ # type: ignore

Partition: TypeAlias = tuple[int, ...]

# Symmetric Functions over Rational Field in the Schur basis
sym_f = SymmetricFunctions(QQ).schur()


class MultiKron(ABC):
    def __init__(self, n: int, N: int):
        pass

    @abstractmethod
    def __call__(self, L: list[Partition]) -> int:
        ...


class MultiKronOrig(MultiKron):
    def __call__(self, L: list[Partition]) -> int :
        """
        L is a list of partitions of the same weight, length at least 2
        return the multiple Kronecker coeffient
        """

        if len(L) == 2:
            if L[1] == L[0]:
                return 1
            else:
                return 0
            
        product = sym_f(list(L[0])).kronecker_product(sym_f(list(L[1])))
        #### This part is unuseful but should go slightly faster with
        if len(L) == 3: # We look for L[2] in product
            for monomial, coeff in product.monomial_coefficients().items():
                if tuple(monomial) == L[2]:
                    return int(coeff)
            return 0
        #### end unuseful
        
        # recursive part
        tot = 0
        for monomial, coeff in product.monomial_coefficients().items():
            tot += coeff * self(L[2:] + [monomial])
        return tot    


class MultiKronOrigFix1(MultiKron):
    """ Enhancement of the original version (using dict.get) """
    def __call__(self, L: list[Partition]) -> int:
        from sage.all import Partition as SagePartition
        assert len(L) >= 2

        if len(L) == 2:
            return 1 if all(lhs == rhs for lhs, rhs in zip(L[0], L[1])) else 0
        
        product = sym_f(tuple(L[0])).kronecker_product(sym_f(tuple(L[1]))).monomial_coefficients()
        if len(L) == 3:
            return int(product.get(SagePartition(L[2]), 0)) # type: ignore
        
        return sum(coeff * self(L[2:] + [monomial]) for monomial, coeff in product.items())


class MultiKronOrigFix2(MultiKron):
    """ Enchancement of the original ersion (using dict.get and intersection of the decomposition) """
    def __call__(self, L: list[Partition]) -> int:
        from sage.all import Partition as SagePartition
        assert len(L) >= 2

        if len(L) == 2:
            return 1 if all(lhs == rhs for lhs, rhs in zip(L[0], L[1])) else 0
        
        product = sym_f(tuple(L[0])).kronecker_product(sym_f(tuple(L[1]))).monomial_coefficients()
        if len(L) == 3:
            return int(product.get(SagePartition(L[2]), 0)) # type: ignore
        
        if len(L) == 4:
            other_product = sym_f(tuple(L[2])).kronecker_product(sym_f(tuple(L[3]))).monomial_coefficients()
            return sum(int(product[m]) * int(other_product[m]) for m in product.keys() & other_product.keys())

        return sum(coeff * self(L[2:] + [monomial]) for monomial, coeff in product.items())


class MultiKronFullProductSage(MultiKron):
    """ Computing the decomposition of the whole product using Sage """
    def __call__(self, L: list[Partition]) -> int:
        from sage.all import Partition as SagePartition
        product = sym_f(tuple(L[0]))
        for p in L[1:-1]:
            product = product.kronecker_product(sym_f(tuple(p)))
        return int(product.monomial_coefficients().get(SagePartition(L[-1]), 0)) # type: ignore


def kronecker_product(partitions: list[Partition]) -> dict[Partition, int]:
    """
    Returns the Kronecker product decomposition (partition, coeff) of two given partitions

    This function needs to be called from Sage (or using sage-python) !
    """
    # Kronecker product
    product = sym_f(tuple(partitions[0]))
    for p in partitions[1:]:
        product = product.kronecker_product(sym_f(tuple(p)))

    # Splitting decomposition
    result: dict[Partition, int] = dict()
    for monomial, coeff in product.monomial_coefficients().items():
        result[tuple(map(int, monomial))] = int(coeff)

    return result


def sort_partitions(partitions: list[Partition]) -> list[Partition]:
    """ Returns sorted nuplets of partitions following lemma 2.2 """
    return sorted(partitions, reverse=True)


class KroneckerCoefficient:
    """ Base class for computing the Kronecker coefficient of a triplet of partitions """

    def product(self, partitions: list[Partition]) -> dict[Partition, int]:
        return kronecker_product(partitions[:-1])
    
    def _kernel(self, partitions: list[Partition]) -> int:
        """ Internal kernel that compute the Kronecker coefficient """
        return self.product(partitions[:-1]).get(partitions[-1], 0)

    def __call__(self, partitions: list[Partition]) -> int:
        partitions = sort_partitions(partitions)
        return self._kernel(partitions)


class KroneckerCoefficientCache(KroneckerCoefficient):
    """ Kronecker coefficient computation using cache to store the computed decompositions """
    def __init__(self, file_prefix: Optional[str] = None):
        super().__init__()
        self.cache = dict()

        if file_prefix is not None:
            self.load_cache(file_prefix)

    def product(self, partitions: list[Partition]) -> dict[Partition, int]:
        try:
            return self.cache[tuple(partitions)]
        except KeyError:
            product = kronecker_product(partitions)
            self.cache[tuple(partitions)] = product
            return product

    def load_cache(self, file_prefix: str) -> None:
        """ Load cache from given filename prefix (will append .pkl.xz) """
        import lzma, pickle, re
        file_name = file_prefix + ".pkl.xz"
        if re.match(r"^https?://", file_name):
            import urllib.request
            with urllib.request.urlopen(file_name) as rfh:
                with lzma.open(rfh, "rb") as fh:
                    self.cache = pickle.load(fh)
        else:
            with lzma.open(file_name, "rb") as fh:
                self.cache = pickle.load(fh)


class MultiKronNupletCache(MultiKron):
    """ Using a cache of the whole nuplet product """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cache = KroneckerCoefficientCache()

    def __call__(self, L: list[Partition]) -> int:
        """ Using a cache of the nuplet of partitions """
        return self.cache(L)


class MultiKronTripletCache(MultiKron):
    """ Using a cache for each triplet """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cache = KroneckerCoefficientCache()

    def __call__(self, L: list[Partition]) -> int:
        """ Using a cache of the nuplet of partitions """
        if len(L) == 2:
            return 1 if all(lhs == rhs for lhs, rhs in zip(L[0], L[1])) else 0
        
        if len(L) == 3:
            return self.cache(L)
        
        if len(L) == 4:
            lhs_product = self.cache.product(L[:2])
            rhs_product = self.cache.product(L[2:])
            return sum(lhs_product[m] * rhs_product[m] for m in lhs_product.keys() & rhs_product.keys())

        product = self.cache.product(L[:2])
        return sum(coeff * self(L[2:] + [monomial]) for monomial, coeff in product.items())    


class MultiKronTripletCacheLoad(MultiKronTripletCache):
    """ Using a cache for each triplet """
    url: str = "https://plmlab.math.cnrs.fr/gt-ia-icj/kronecker/-/raw/main/data/KroneckerCoefficientCache_Sage_n{}.pkl.xz?ref_type=heads&inline=false"
    def __init__(self, n: int, N: int):
        super().__init__(n, N)
        self.cache.load_cache(self.url.format(n))


def partition_generator(N: int) -> Iterable[Partition]:
    """
    Generates all partitions of an integer N >= 0.

    Could be optimized but it is clearly enough for the N we will consider.
    """
    if N <= 0:
        yield ()
        return

    for head in range(N, 0, -1):
        for tail in partition_generator(N - head):
            if len(tail) == 0 or head >= tail[0]:
                yield (head,) + tail


class AllNupletOfFixedN(Iterable):
    """ Returns an iterable over all N-uplets of partitions of n

    Supports len function.
    """
    partitions: list[Partition]
    length: int
    N: int

    def __init__(self, n: int, N: int):
        self.partitions = list(partition_generator(n))
        self.length = len(self.partitions)**N
        self.N = N

    def __len__(self) -> int:
        return self.length

    def __iter__(self) -> Iterator[list[Partition]]:
        return (list(nuplet) for nuplet in itertools.product(self.partitions, repeat=self.N))
    

if __name__ == "__main__":
    all_methods = ["Orig", "OrigFix1", "OrigFix2", "FullProductSage", "NupletCache", "TripletCache", "TripletCacheLoad"]

    import argparse
    parser = argparse.ArgumentParser(
        "Comparing methods to compute Kronecker coefficients of multiple partitions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("n", type=int, help="The length of the partitions")
    parser.add_argument("N", type=int, default=3, help="Computing the Kronecker coefficient for a N-uplet of partitions")
    parser.add_argument("--methods", choices=all_methods, nargs="*", default=all_methods, help="Methods to compare")
    config = parser.parse_args()

    print(f"Computing Kronecker coefficient of {config.N}-uplet of partitions of {config.n}:")

    all_nuplets = list(AllNupletOfFixedN(config.n, config.N))

    import random
    partitions_check = random.sample(all_nuplets, 10)
    partitions_check = []
    print("Reference partitions:")
    for p in partitions_check:
        print(f"\t{p}")

    from tqdm.auto import tqdm
    for method in config.methods:
        print(f"\nUsing method {method}:")
        fn = globals()[f"MultiKron{method}"](config.n, config.N)

        for nuplet in tqdm(all_nuplets):
            fn(nuplet)

        print("Check on the reference partitions:", [fn(p) for p in partitions_check])
