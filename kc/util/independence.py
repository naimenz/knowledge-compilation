"""File for the independence compilation rule"""

from kc.data_structures import *
from kc.util import *

def independence(delta: 'CNF') -> 'CNF':
    """Apply the independence algorithm to a theory delta
    (plus check if it is possible).

    Returns an extensional conjunction of two independent subtheories"""

    subtheory_pairs = generate_all_subtheory_pairs(delta)
    for pair in subtheory_pairs:
        delta1, delta2 = pair
        if theories_independent(delta1, delta2):
            return compile_theory(delta1).join(compile_theory(delta2))

