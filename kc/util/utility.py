"""This is a file for *small* utility functions
that are used in various places"""

# for powerset generation
from itertools import chain, combinations

from typing import TypeVar, Iterable, Set, FrozenSet, Tuple, Sequence, Generator, List, Any, Union
T = TypeVar('T')

def powerset(iterable: Iterable[T]) -> Iterable[Tuple[T, ...]]:
    """Modified from https://stackoverflow.com/a/1482316/10005793
    Returns all non-empty subsets.
    e.g. powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def get_element_of_set(s: Iterable[T]) -> T:
        """Simple function to get a (random) element of a set without mutating it.
        Useful for getting the only element of a 1-element set."""
        return next(iter(s))


def partition_set(s: Set[T]) -> Generator[Set[FrozenSet[T]], None, None]:
    """Wrapper around the partition function to return a set of sets"""
    collection = list(s)
    partitioned_lists = _partition(collection)
    for partition_list in partitioned_lists:
        yield set(frozenset(li) for li in partition_list)


def _partition(collection: List[T]) -> Generator[List[List[T]], None, None]:
    """Function to partition a list into lists of lists.
    From https://stackoverflow.com/a/30134039/10005793."""
    if len(collection) <= 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in _partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller
