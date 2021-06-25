"""File for vacuous conjunction compilation rule"""

from kc.data_structures import ForAllNode
from kc.compiler import KCRule

from typing import Tuple
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.data_structures import CNF
    from kc.compiler import Compiler

class VacuousConjunction(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, None]:
        """VacuousConjunction is applicable if the theory consists of
        a single clause with no bound variables.
        Returns a boolean plus None, because no stored data is needed"""
        applicable = len(cnf.clauses) == 1 and len(tuple(cnf.clauses)[0].bound_vars) == 0
        return applicable, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: None, compiler: 'Compiler') -> 'ForAllNode':
        """Apply VacuousConjunction and return an NNFNode
        TODO: Decide whether accepting UnconstrainedClauses in compile
        is the right approach."""
        clause = tuple(CNF.clauses)[0] # only one so we can access it directly
        child = compiler.compile(clause.unconstrained_clause)
        return ForAllNode(child, clause.bound_vars, clause.cs)


