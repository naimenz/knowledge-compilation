"""File for vacuous conjunction compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple

class VacuousConjunction(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, None]:
        """VacuousConjunction is applicable if the theory consists of
        a single clause with no bound variables.
        Returns a boolean plus None, because no stored data is needed"""
        applicable = len(cnf.clauses) == 1 and len(list(cnf.clauses)[0].bound_vars) == 0
        return applicable, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: None) -> 'NNFNode':
        """Apply VacuousConjunction and return an NNFNode"""
        raise NotImplementedError('VacuousConjunction.apply not implemented')

