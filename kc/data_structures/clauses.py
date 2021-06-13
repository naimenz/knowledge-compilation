"""
Classes for clauses in FOL-DC.
This includes constrained AND unconstrained clauses.

TODO: Figure out if the inheritance structure for UnitClause and ConstrainedAtom makes sense.
"""

from kc.data_structures.logicalterms import *
from kc.data_structures.literals import *
from kc.data_structures.constraints import *
from abc import ABC

from typing import List

class Clause(ABC):
    """Abstract base class for constrained and unconstrained clauses"""

class UnconstrainedClause(Clause):
    """An FOL unconstrained clause.
    This consists of a set of FOL literals, which form a disjunction
    """

    def __init__(self, literals: Iterable['Literal']) -> None:
        """NOTE: using a set to represent literals.
        This is justified because the order is unimportant and repeated literals 
        in a clause are redundant (they cannot change the disjunction)"""
        self.literals = frozenset(literals)

    def __eq__(self, other: Any) -> bool:
        """Two unconstrained clauses are equal if they have the same literals

        NOTE: for now the ordering has to be the same
        TODO: make literals hashable so I can compare as sets"""
        if not isinstance(other, UnconstrainedClause):
            return False
        if len(self.literals) != len(other.literals): 
            return False
        same_literals = all(self_l == other_l for self_l, other_l in zip(self.literals, other.literals))
        return same_literals

    def __hash__(self) -> int:
        return hash(self.literals)

    def __str__(self) -> str:
        literal_strs = [str(literal) for literal in self.literals]
        logical_or_string = ' \u2228 '
        return f"({logical_or_string.join(literal_strs)})"

    def __repr__(self) -> str:
        return self.__str__()

class ConstrainedClause(Clause):
    """An FOL-DC constrained clause.
    This consists of an unconstrained clause with a set of bound variables and a constraint set.

    NOTE: For now we only work with bound *logical* variables."""

    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        self.unconstrained_clause = unconstrained_clause
        self.bound_vars = frozenset(bound_vars)
        self.cs = cs

    def __eq__(self, other: Any) -> bool:
        """Two constrained literals are equal if they have the same unconstrained literals, the same constraint sets,
         and the same bound variables"""
        if not isinstance(other, ConstrainedClause):
            return False
        same_u_clause = (self.unconstrained_clause == other.unconstrained_clause)
        if len(self.bound_vars) != len(other.bound_vars): 
            return False
        same_bound_vars = all(self_v == other_v for self_v, other_v in zip(self.bound_vars, other.bound_vars))
        same_cs = (self.cs == other.cs)
        return same_u_clause and same_bound_vars and same_cs

    def __hash__(self) -> int:
       return hash((self.unconstrained_clause, self.bound_vars, self.cs))

    def __str__(self) -> str:
        bound_vars_strs = [str(var) for var in self.bound_vars]
        for_all_string = '\u2200'
        return f"{for_all_string}{{{', '.join(bound_vars_strs)}}}, {self.cs} : {self.unconstrained_clause}"

    def __repr__(self) -> str:
        return self.__str__()


class UnitClause(ConstrainedClause):
    """AN FOL-DC unit clause.
    This is a constrained clause with a single literal
    """
    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(unconstrained_clause.literals) == 1) # ensure that this is a unit clause
        self.literal = list(unconstrained_clause.literals)[0] # convert to 1-item list and get item
        super(UnitClause, self).__init__(unconstrained_clause, bound_vars, cs)


class ConstrainedAtom(UnitClause):
    """An FOL-DC constrained atom.
    This is a constrained clause with a single atom (positive literal).
    """
    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(unconstrained_clause.literals) == 1) # ensure that this is a unit clause
        super(ConstrainedAtom, self).__init__(unconstrained_clause, bound_vars, cs)
        assert(self.literal.polarity) # ensure that it is not negated

    @property
    def atom(self) -> 'Atom':
        """Get the atom from the unconstrained clause"""
        return self.literal.atom


