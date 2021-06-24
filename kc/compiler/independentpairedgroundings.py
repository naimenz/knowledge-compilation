"""File for independent Paired groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class IndependentPairedGroundings(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['EquivalenceClass']]:
        """IndependentPairedGroundings is applicable if the theory is shattered
        (which is represented by a flag) and there is a root unifying class with two variables per clause.
        Returns True and the root unifying class if applicable, and False, None otherwise."""
        needs_shattering = not cnf.shattered
        unifying_classes = cnf.get_unifying_classes()
        # writing a little function to check if root in this particular cnf
        is_root_in_cnf = lambda eq_class: eq_class.is_root_in_cnf(cnf)
        root_unifying_classes = filter(is_root_in_cnf, unifying_classes)
        for root_unifying_class in root_unifying_classes:
            if cnf.eq_class_has_two_variables(root_unifying_class):
                return True, root_unifying_class
        return False, None


    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: 'EquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentPairedGroundings and return an NNFNode"""
        raise NotImplementedError('IndependentPairedGroundings.apply not implemented')

