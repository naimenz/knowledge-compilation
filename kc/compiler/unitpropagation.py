"""File for unit propagation compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Optional, Tuple
class UnitPropagation(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional['ConstrainedClause']]:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)
        Returns True and the unit clause if applicable, and False, None otherwise"""
        for clause in delta.clauses:
            if len(clause.unconstrained_clause.literals) == 1:
                return True, clause
        return False, None

    @classmethod
    def apply(cls, delta: 'CNF', stored_data: 'ConstrainedClause') -> 'NNFNode':
        """Apply UnitPropagation (which requires applying splitting and conditioning)
        and return an NNFNode"""
        raise NotImplementedError('UnitPropagation.apply not implemented')

