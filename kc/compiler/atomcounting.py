"""File for atom counting compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_constrained_atoms

from typing import Tuple, Optional

class AtomCounting(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional['ConstrainedAtom']]:
        """AtomCounting is applicable if the theory is shattered
        (which is represented by a flag) and there is an atom in delta with exactly one bound
        logical variable
        Returns True and the atom if applicable, and False, None otherwise."""
        needs_shattering = not delta.shattered
        if needs_shattering:
            return False, None
        for clause in delta.clauses:
            for c_atom in get_constrained_atoms(clause):
                overlap = set(c_atom.atom.terms).intersection(c_atom.bound_vars) 
                if len(overlap) == 1:
                    return True, c_atom
        return False, None

    @classmethod
    def apply(cls, delta: 'CNF', stored_data: 'ConstrainedAtom') -> 'NNFNode':
        """Apply AtomCounting and return an NNFNode"""
        raise NotImplementedError('AtomCounting.apply not implemented')

