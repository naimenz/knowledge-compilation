"""File for atom counting compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class AtomCounting(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['ConstrainedAtom']]:
        """AtomCounting is applicable if the theory is shattered
        (which should have been checked before this rule) and there is an atom in cnf with exactly one bound
        logical variable
        Returns True and the atom if applicable, and False, None otherwise."""
        # TODO: Heuristic for deciding which c_atom to use (from Forclift)
        for clause in cnf.clauses:
            for c_atom in clause.get_constrained_atoms():
                overlap = set(c_atom.atom.terms).intersection(c_atom.bound_vars) 
                if len(overlap) == 1:
                    return True, c_atom
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: 'ConstrainedAtom', compiler: 'Compiler') -> 'NNFNode':
        """Apply AtomCounting and return an NNFNode"""
        raise NotImplementedError('AtomCounting.apply not implemented')

