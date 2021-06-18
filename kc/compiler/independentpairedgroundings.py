"""File for independent Paired groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple

class IndependentPairedGroundings(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, None]:
        """IndependentPairedGroundings is applicable if the theory is shattered
        (which is represented by a flag) and there is a root unifying class with two variables per clause.
        Returns True and the root unifying class if applicable, and False, None otherwise."""
        needs_shattering = not delta.shattered

    @classmethod
    def apply(cls, delta: 'CNF', stored_data: None) -> 'NNFNode':
        """Apply IndependentPairedGroundings and return an NNFNode"""
        raise NotImplementedError('IndependentPairedGroundings.apply not implemented')

