"""
Class for FO-CNF formulas.
"""

from kc.data_structures import ConstrainedClause, EquivalenceClass

from typing import List, Any, Iterable

class CNF:
    """
    A FOL-DC CNF.
    This consists of a set of constrained clauses, which form a conjunction.
    """

    def __init__(self, clauses: Iterable['ConstrainedClause']) -> None:
        self.clauses = frozenset(clauses)
        self.shattered = False # keep track of whether this cnf has undergone shattering

    def join(self, other: 'CNF') -> 'CNF':
        """Combine two CNFs into one.
        TODO: Eventually this will probably be done with circuit nodes and/or sets"""
        return CNF(self.clauses.union(other.clauses))

    # def apply_substitution(self, substitution: 'Substitution') -> 'CNF':
    #     """Return a new CNF, the result of applying substitution to this CNF"""
    #     new_clauses = set(clause.apply_substitution(substitution) for clause in self.clauses)
    #     return CNF(new_clauses)

    def eq_class_has_one_variable(self, eq_class: 'EquivalenceClass') -> bool:
        """Determine whether a given root equivalence class has a single bound variable
        per clause or not.
        NOTE: We assume that this equivalence class is root in the cnf"""
        for clause in self.clauses: 
            if len(eq_class.members.intersection(clause.bound_vars)) != 1:
                return False
        return True

    def eq_class_has_two_variables(self, eq_class: 'EquivalenceClass') -> bool:
        """Determine whether a given root equivalence class has two bound variables
        per clause or not.
        NOTE: We assume that this equivalence class is root in the cnf"""
        for clause in self.clauses: 
            if len(eq_class.members.intersection(clause.bound_vars)) != 2:
                return False
        return True

    def __eq__(self, other: Any) -> bool:
        """Two CNFs are equal if they have the same clauses"""
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

