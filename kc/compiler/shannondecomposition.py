"""File for shannon decomposition compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler
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
    def apply(cls, cnf: 'CNF', unbound_atom: 'ConstrainedAtom', compiler: 'Compiler') -> 'NNFNode':
        """Apply ShannonDecomposition and return an NNFNode"""
        true_literal = unbound_atom.to_unit_clause()
        false_literal = UnitClause([Literal(unbound_atom.atom, False)], unbound_atom.bound_vars, unbound_atom.cs)
        true_branch = cnf.join(CNF([true_literal]))
        false_branch = cnf.join(CNF([false_literal]))

        return AndNode(compiler.compile(true_branch), compiler.compile(false_branch))

