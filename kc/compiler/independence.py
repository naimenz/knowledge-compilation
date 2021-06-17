"""File for independence compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import partition_cnf

class Independence(KCRule):
    
    @staticmethod
    def is_applicable(delta: 'CNF') -> bool:
        """Independence is applicable if the theory can be divided into
        two subtheories such that the substheories make up the whole theory and are
        independent."""
        for delta1, delta2 in partition_cnf(delta):
            if are_independent_clauses(delta1, delta2):
                return True
        return False


    @staticmethod
    def apply(delta: 'CNF') -> 'NNFNode':
        """Apply VacuousConjunction and return an NNFNode"""
        raise NotImplementedError('VacuousConjunction.apply not implemented')

