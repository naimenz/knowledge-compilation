"""File for leaf construction. This isn't formally a compilation rule from the PhD,
but it is described in the Compile algorithm there. For consistency, we make it into a 
full rule."""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class CheckTautology(KCRule):
    """To simplify situations where one clause in a theory is a tautology,
    we check whether each clause is a tautology BEFORE AtomCounting (in an attempt to fix
    a bug where AtomCounting is done unnecessarily)."""
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['Clause']]:
        for clause in sorted(cnf.clauses):
            if clause.is_tautology():
                return True, clause
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', tautology_clause: 'Clause', compiler: 'Compiler') -> 'NNFNode':
        """Remove tautology and create new CNF without it."""
        simplified_cnf = CNF(cnf.clauses.difference({tautology_clause}),
                             shattered=cnf.shattered,
                             subdivided=cnf.subdivided)
        return AndNode(compiler.compile(simplified_cnf), TrueNode())
