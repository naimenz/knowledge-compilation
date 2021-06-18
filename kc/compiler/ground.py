"""File for ground compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_logical_variables_from_cs

from typing import Tuple, Optional

class Ground(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional['ConstraintSet']]:
        """Ground is applicable if the theory is contains a constraint set with no free variables.
        Returns True and the constraint set if applicable, and False, None otherwise.
        NOTE: for now, this only considers LogicalVariables, not DomainVariables"""
        for clause in delta.clauses:
            cs_variables = get_logical_variables_from_cs(clause.cs)
            variable_overlap = clause.bound_vars.intersection(cs_variables)
            if variable_overlap == cs_variables: # all variables in cs are in bound
                return True, clause.cs
        return False, None

    @classmethod
    def apply(cls, delta: 'CNF', stored_data: 'ConstraintSet') -> 'NNFNode':
        """Apply Ground and return an NNFNode"""
        raise NotImplementedError('Ground.apply not implemented')

