"""File for shatter compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple

class ShatteredCompilation(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, None]:
        """ShatteredCompilation is applicable if the theory is not already shattered (which is represented by a flag )
        Returns True or False depending on the flag, plus None (no stored data needed)"""
        shattering_applicable = not cnf.shattered
        return shattering_applicable, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: None) -> 'NNFNode':
        """Apply ShatteredCompilation and return an NNFNode"""
        raise NotImplementedError('ShatteredCompilation.apply not implemented')

