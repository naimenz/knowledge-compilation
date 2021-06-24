"""File for independent single groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class IndependentSingleGroundings(KCRule):
    
    # TODO: Try checking ISG and IPG together
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['EquivalenceClass']]:
        """IndependentSingleGroundings is applicable if the theory is shattered
        (which it must already be to each this rule) and there is a root unifying class with one variable per clause.
        Returns True and the root unifying class if applicable, and False, None otherwise."""
        unifying_classes = cnf.get_unifying_classes()
        # writing a little function to check if root in this particular cnf
        # TODO: Think about making this more readable or faster
        is_root_in_cnf = lambda eq_class: eq_class.is_root_in_cnf(cnf)
        root_unifying_classes = filter(is_root_in_cnf, unifying_classes)
        for root_unifying_class in root_unifying_classes:
            if cnf.eq_class_has_one_variable(root_unifying_class):
                return True, root_unifying_class
        return False, None


    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: 'EquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentSingleGroundings and return an NNFNode"""
        raise NotImplementedError('IndependentSingleGroundings.apply not implemented')

