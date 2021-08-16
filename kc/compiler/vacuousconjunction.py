"""File for vacuous conjunction compilation rule"""

from kc.data_structures import ForAllNode, UnconstrainedClause, CNF
from kc.compiler import KCRule

from typing import Tuple
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class VacuousConjunction(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, None]:
        """VacuousConjunction is applicable if the theory consists of
        a single CONSTRAINED clause with no bound variables (i.e. it must still have a constraint set).
        Returns a boolean plus None, because no stored data is needed"""
        if not (len(cnf.c_clauses) == 1 and len(cnf.u_clauses) == 0):
            return False, None

        if len(tuple(cnf.c_clauses)[0].bound_vars) == 0:
            return True, None
        else:
            return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: None, compiler: 'Compiler') -> 'ForAllNode':
        """Apply VacuousConjunction and return an NNFNode"""
        clause = tuple(cnf.c_clauses)[0] # only one so we can access it directly
        child_cnf = CNF([UnconstrainedClause(clause.literals)], subdivided=cnf.subdivided)
        child = compiler.compile(child_cnf)
        return ForAllNode(child, clause.bound_vars, clause.cs)
