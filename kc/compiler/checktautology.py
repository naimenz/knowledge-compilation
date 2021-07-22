"""File for leaf construction. This isn't formally a compilation rule from the PhD,
but it is described in the Compile algorithm there. For consistency, we make it into a 
full rule."""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class CheckTautology(KCRule):
    """To simplify situations where one clause in a theory is a tautology,
    we check each whether each clause is a tautology BEFORE AtomCounting (in an attempt to fix
    a bug where AtomCounting is done unnecessarily)."""
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['Clause']]:
        """LeafConstruction is applicable if the theory is a single True, False, or literal.
        Returns True or False depending on the flag, plus None (no stored data needed)
        NOTE: I'm not quite sure what it means for a clause to be simply True
        TODO: Figure out the logic here instead of just returning False"""
        for clause in sorted(cnf.clauses):
            if clause.is_tautology():
                return True, clause
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', tautology_clause: 'Clause', compiler: 'Compiler') -> 'NNFNode':
        """Apply LeafConstruction and return an NNFNode"""
        simplified_cnf = CNF(cnf.clauses.difference(set([tautology_clause])), shattered=cnf.shattered)
        return AndNode(compiler.compile(simplified_cnf), TrueNode())
