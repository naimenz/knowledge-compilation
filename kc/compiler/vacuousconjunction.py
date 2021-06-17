"""File for vacuous conjunction compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

class VacuousConjunction(KCRule):
    
    @staticmethod
    def is_applicable(delta: 'CNF') -> bool:
        """VacuousConjunction is applicable if the theory consists of
        a single clause with no bound variables."""
        return len(delta.clauses) == 1 and len(list(delta.clauses)[0].bound_vars) == 0

    @staticmethod
    def apply(delta: 'CNF') -> 'NNFNode':
        """Apply VacuousConjunction and return an NNFNode"""
        raise NotImplementedError('VacuousConjunction.apply not implemented')

