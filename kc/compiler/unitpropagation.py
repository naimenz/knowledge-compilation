"""File for unit propagation compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

class UnitPropagation(KCRule):
    
    @staticmethod
    def is_applicable(delta: 'CNF') -> bool:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)"""
        for clause in delta.clauses:
            if len(clause.unconstrained_clause.literals) == 1:
                return True
        return False

    @staticmethod
    def apply(delta: 'CNF') -> 'NNFNode':
        """Apply UnitPropagation (which requires applying splitting and conditioning)
        and return an NNFNode"""
        raise NotImplementedError('UnitPropagation.apply not implemented')

