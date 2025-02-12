# Goal

The goal is to compare different partitions sorting strategy when computing and caching decompositions of Kronecker product of multiple partitions.

# Implementation

The reference implementation of the computation of the Kronecker product is the class `KroneckerCoefficientMLCache`
that relies on a divide-to-conquer strategy on the list of partitions so that to reduce the number
of products in cache.

# Problems and test codes
Two problems are considered:

1. computing the Kronecker coefficient of the 11390625 6-uplets of partitions of 7 using the code
    ```Python
    kc = KroneckerCoefficientMLCache()
    all_nuplets = (tuple(nuplet) for nuplet in itertools.product(Partition.all_for_integer(7), repeat=6))
    for nuplet in tqdm(all_nuplets):
        kc(nuplet)
    ```
2. generating all 4-uplets of partitions of 11 that have a non-zero Kronecker coefficient and with constraints on the length of each partition. This is tested using the code
    ```Python
    kc = KroneckerCoefficientMLCache()
    ascending = all_partitions_of_max_length(11, [2, 3, 4, 5], kc)
    descending = all_partitions_of_max_length(11, [5, 4, 3, 2], kc)
    for _ in tqdm(itertools.chain(ascending, descending)): pass
    ```

After launching the test code, we consider multiples metrics :

1. the computation time in seconds displayed by `tqdm`,
2. the total number of products in the cache,
3. the number of products in the cache for each number of partitions multiplied together. 

These two last informations are gathered by displaying the variable `kc`.

# Sorting strategies

Multiple sorting strategies are compared, at two different levels.

In the Kronecker computation class, the N-uplet of partitions is sorted before
checking the cache with the $N-1$ first elements using the following strategies:

1. no sorting at all,
2. sorting by lexicographical order,
3. sorting by length (increasing and decreasing) only,
4. sorting by length and lexicographical order (increasing and decreasing).

When generating all non-zero nuplets with length constraints, we check the following strategies of sorting the constraints :

1. no sorting at all,
2. ascending order,
3. descending order.

# Results

| nuplet sort | constraint sort | all time | all cache | constraints time | constraints cache |
| --- | --- | --- | --- | --- | --- |
| none | none | 105s | 3600 (#2=225,#3=3375) | 665s | 18039 (#2=999,#3=17040) |
| none | asc  | - | - | 73s | 2688 (#2=96,#3=2592) |
| none | desc | - | - | 637s | 16983 (#2=999,#3=15984) |
| lexico asc | none | 97s | 800 (#2=120,#3=680) | 218s | 6968 (#2=480,#3=6488) |
| lexico asc | asc | - | - | 40s | 1628 (#2=141,#3=1487) |
| lexico asc | desc | - | - | 222s | 6968 (#2=480,#3=6488) |
| lexico desc | none | 99s | 800 (#2=120,#3=680) | 306s | 7136 (#2=648,#3=6488)
| lexico desc | asc | - | - | 65s | 1799 (#2=312,#3=1487) |
| lexico desc | desc | - | - | 311s | 7136 (#2=648,#3=6488)
| length asc | none | 100s | 1046 (#2=133,#3=913) | 340s | 10140 (#2=367,#3=9773) |
| length asc | asc | - | - | 53s | 2158 (#2=91,#3=2067) |
| length asc | desc | - | - | 336s | 10140 (#2=367,#3=9773) |
| length desc | none | 96s | 1046 (#2=133,#3=913) | 436s | 10531 (#2=758,#3=9773) |
| length desc | asc | - | - | 85s | 2434 (#2=367,#3=2067)
| length desc | desc | - | - | 443s | 10531 (#2=758,#3=9773) |
| len & lexico asc | none | **100s** | **800 (#2=120,#3=680)** | 219s | 6800 (#2=312,#3=6488) |
| **len & lexico asc** | **asc** | - | - | **39s** | **1568 (#2=81,#3=1487)** |
| len & lexico asc | desc | - | - | 222s | 6800 (#2=312,#3=6488) |
| len & lexico desc | none | 102s | 800 (#2=120,#3=680) | 305s | 7136 (#2=648,#3=6488) |
| len & lexico desc | asc | - | - | 64s | 1799 (#2=312,#3=1487) |
| len & lexico desc | desc | - | - | 305s | 7136 (#2=648,#3=6488) |
