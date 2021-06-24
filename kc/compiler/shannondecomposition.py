"""File for shannon decomposition compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional

class ShannonDecomposition(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['ConstrainedAtom']]:
        """ShannonDecomposition is applicable if the theory contains
        an atom without bound logical variables.
        Returns True plus the atom if applicable, and False, None otherwise"""
        for clause in cnf.clauses:
            for c_atom in clause.get_constrained_atoms():
                overlap = c_atom.bound_vars.intersection(set(c_atom.atom.terms)) 
                if len(overlap) == 0:
                    return True, c_atom
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: 'ConstrainedAtom') -> 'NNFNode':
        """Apply ShannonDecomposition and return an NNFNode"""
        raise NotImplementedError('ShannonDecomposition.apply not implemented')

