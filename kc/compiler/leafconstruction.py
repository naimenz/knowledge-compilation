"""File for leaf construction. This isn't formally a compilation rule from the PhD,
but it is described in the Compile algorithm there. For consistency, we make it into a 
full rule."""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple

class LeafConstruction(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, None]:
        """LeafConstruction is applicable if the theory is a single True, False, or literal.
        Returns True or False depending on the flag, plus None (no stored data needed)
        NOTE: I'm not quite sure what it means for a clause to be simply True
        TODO: Figure out the logic here instead of just returning False"""
        return False, None

    @classmethod
    def apply(cls, delta: 'CNF', stored_data: None) -> 'NNFNode':
        """Apply LeafConstruction and return an NNFNode"""
        raise NotImplementedError('LeafConstruction.apply not implemented')

