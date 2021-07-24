
from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from copy import deepcopy

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

def smooth(root: 'NNFNode', cnf: 'CNF') -> 'NNFNode':
    """This is not a compilation rule so much as a post-processing step.
    We pass in the input CNF so that we know which predicates and domains need to be covered.
    It is applied recursively, starting with the root NNFNode and smoothing down each branch.
    It returns a smoothed version of the NNF below the input node
    """
    starting_node = deepcopy(root)
