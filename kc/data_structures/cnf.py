"""
Class for FO-CNF formulas.
"""

from kc.data_structures.clauses import *

from typing import List, Any

class FO_CNF:
    """
    A FOL-DC CNF.
    This consists of a set of constrained clauses, which form a conjunction.
    """

    def __init__(self, clauses: List['ConstrainedClause']) -> None:
        self.clauses = clauses

    def __eq__(self, other: Any) -> bool:
        """Two CNFs are equal if they have the same clauses

        NOTE: for now the ordering of the clauses is important.
        TODO: make clauses hashable so I can compare sets"""
        if not isinstance(other, FO_CNF):
            return False
        same_clauses = all(self_clause == other_clause for self_clause, other_clause in zip(self.clauses, other.clauses))
        return same_clauses

    def __str__(self) -> str:
        clause_strs = [f'({str(clause)})' for clause in self.clauses]
        return '\nAND\n'.join(clause_strs)

    def __repr__(self) -> str:
        return self.__str__()

