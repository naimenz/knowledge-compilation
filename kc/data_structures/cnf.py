"""
Class for FO-CNF formulas.
"""

from kc.data_structures.clauses import *

from typing import List, Any

class CNF:
    """
    A FOL-DC CNF.
    This consists of a set of constrained clauses, which form a conjunction.
    """

    def __init__(self, clauses: Iterable['ConstrainedClause']) -> None:
        """TODO: replace list with set"""
        self.clauses = frozenset(clauses)
        self.shattered = False # keep track of whether this cnf has undergone shattering

    def join(self, other: 'CNF') -> 'CNF':
        """Combine two CNFs into one.
        TODO: Eventually this will probably be done with circuit nodes and/or sets"""
        return CNF(self.clauses.union(other.clauses))

    def apply_substitution(self, substitution: 'Substitution') -> 'CNF':
        """Return a new CNF, the result of applying substitution to this CNF"""
        new_clauses = set(clause.apply_substitution(substitution) for clause in self.clauses)
        return CNF(new_clauses)

    def __eq__(self, other: Any) -> bool:
        """Two CNFs are equal if they have the same clauses

        NOTE: for now the ordering of the clauses is important.
        TODO: make clauses hashable so I can compare sets"""
        if not isinstance(other, CNF):
            return False
        same_clauses = (self.clauses == other.clauses)
        return same_clauses

    def __hash__(self) -> int:
        return hash(self.clauses)

    def __str__(self) -> str:
        clause_strs = [f'({str(clause)})' for clause in self.clauses]
        return '\nAND\n'.join(clause_strs)

    def __repr__(self) -> str:
        return self.__str__()

