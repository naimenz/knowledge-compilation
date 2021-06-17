"""This file contains a function to partition a set into k-sized blocks.
The code is from stackoverflow and is a bit difficult to read.
Link: https://codereview.stackexchange.com/a/240277

I've written a small wrapper around it to make it work with CNFs directly
"""
from kc.data_structures import *

from typing import Tuple, Generator

def partition_k(collection, minimum, k):
  if len(collection) == 1:
    yield [ collection ]
    return

  first = collection[0]
  for smaller in partition_k(collection[1:], minimum - 1, k):
    if len(smaller) > k: continue
    # insert `first` in each of the subpartition's subsets
    if len(smaller) >= minimum:
      for n, subset in enumerate(smaller):
        yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
    # put `first` in its own subset 
    if len(smaller) < k: yield [ [ first ] ] + smaller


def partition_cnf(cnf: 'CNF') -> Generator[Tuple['CNF', 'CNF'], None, None]:
    """Return a generator that produces all possible splits of a CNF into
    two non-empty sets of clauses"""
    all_clauses = list(cnf.clauses)
    partition_generator = partition_k(all_clauses, 2, 2)
    for partition in partition_generator:
        yield (CNF(partition[0]), CNF(partition[1]))
    return
