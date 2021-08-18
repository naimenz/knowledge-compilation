"""File for ground compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class Ground(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['ConstrainedClause']]:
        """Ground is applicable if the theory is contains a constraint set with no free variables.
        Returns True and the clause with that constraint set if applicable, and False, None otherwise.
        NOTE: for now, this only considers LogicalVariables, not DomainVariables"""
        for clause in cnf.c_clauses:
            cs_variables = clause.cs.get_logical_variables()
            variable_overlap = clause.bound_vars.intersection(cs_variables)
            if variable_overlap == cs_variables:  # all variables in cs are in bound
                return True, clause.cs
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', clause_to_ground: 'ConstrainedClause', compiler: 'Compiler') -> 'NNFNode':
        """Apply Ground and return an NNFNode"""
        raise NotImplementedError('Ground.apply not implemented')

