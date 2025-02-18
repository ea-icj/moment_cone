__all__ = (
    "unique_combinations",
)

def unique_combinations(mylist, k):
    """
    The list is sorted.
    The list is viewed as a multiset. 
    Create the list of multisets contained in list of cardinality k

    Example :
    >>> unique_combinations([3,3,3,2,2,1],3)
    [(3, 3, 3), (3, 3, 2), (3, 3, 1), (3, 2, 2), (3, 2, 1), (2, 2, 1)]
    """
    def backtrack(start, current):
        if len(current) == k:
            result.append(tuple(current))
            return
        for i in range(start, len(mylist)):
            # Éviter de répéter un élément au même niveau
            if i > start and mylist[i] == mylist[i - 1]:
                continue
            current.append(mylist[i])
            backtrack(i + 1, current)
            current.pop()

    result = []
    backtrack(0, [])
    return result

