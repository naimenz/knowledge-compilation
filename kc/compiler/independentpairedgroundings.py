"""File for independent Paired groundings compilation rule"""

from kc.data_structures import *
from kc.util import *
from kc.compiler import KCRule

from typing import Tuple, Optional

class IndependentPairedGroundings(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional['EquivalenceClass']]:
        """IndependentPairedGroundings is applicable if the theory is shattered
        (which is represented by a flag) and there is a root unifying class with two variables per clause.
        Returns True and the root unifying class if applicable, and False, None otherwise."""
        needs_shattering = not delta.shattered
        unifying_classes = get_unifying_classes(delta)
        # writing a little function to check if root in this particular cnf
        is_root_in_delta = lambda eq_class: is_root_eq_class(eq_class, delta)
        root_unifying_classes = filter(is_root_in_delta, unifying_classes)
        for root_unifying_class in root_unifying_classes:
            if eq_class_has_two_variables(root_unifying_class, delta):
                return True, root_unifying_class
        return False, None


    @classmethod
    def apply(cls, delta: 'CNF', stored_data: 'EquivalenceClass') -> 'NNFNode':
        """Apply IndependentPairedGroundings and return an NNFNode"""
        raise NotImplementedError('IndependentPairedGroundings.apply not implemented')

