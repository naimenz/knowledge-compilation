"""File for shannon decomposition compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_constrained_atoms

from typing import Tuple, Optional

class ShannonDecomposition(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional['ConstrainedAtom']]:
        """ShannonDecomposition is applicable if the theory contains
        an atom without bound logical variables.
        Returns True plus the atom if applicable, and False, None otherwise"""
        for clause in delta.clauses:
            for c_atom in get_constrained_atoms(clause):
                overlap = c_atom.bound_vars.intersection(set(c_atom.atom.terms)) 
                if len(overlap) == 0:
                    return True, c_atom
        return False, None

    @classmethod
    def apply(cls, delta: 'CNF', stored_data: 'ConstrainedAtom') -> 'NNFNode':
        """Apply ShannonDecomposition and return an NNFNode"""
        raise NotImplementedError('ShannonDecomposition.apply not implemented')

