"""This is a file for *small* utility functions
that are used in various places"""

# for powerset generation
from itertools import chain, combinations

from typing import TypeVar, Iterable, Set, Tuple
T = TypeVar('T')

def powerset(iterable: Iterable[T]) -> Iterable[Tuple[T, ...]]:
    """Modified from https://stackoverflow.com/a/1482316/10005793
    Returns all non-empty subsets.
    e.g. powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def get(s: Set[T]) -> T:
        """Simple function to get a (random) element of a set without mutating it.
        Useful for getting the only element of a 1-element set."""
        return next(iter(s))

            



