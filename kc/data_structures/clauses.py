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

    def __init__(self, literals: List['Literal']) -> None:
        self.literals = literals

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
            bound_vars: List['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        self.unconstrained_clause = unconstrained_clause
        self.bound_vars = bound_vars
        self.cs = cs

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
            bound_vars: List['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(unconstrained_clause.literals) == 1) # ensure that this is a unit clause
        super(UnitClause, self).__init__(unconstrained_clause, bound_vars, cs)


class ConstrainedAtom(UnitClause):
    """An FOL-DC constrained atom.
    This is a constrained clause with a single atom (positive literal).
    """
    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: List['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(unconstrained_clause.literals) == 1) # ensure that this is a unit clause
        assert(unconstrained_clause.literals[0].polarity) # ensure that it is not negated
        super(ConstrainedAtom, self).__init__(unconstrained_clause, bound_vars, cs)


if __name__ == '__main__':
    pred1 = Predicate('smokes', 1)
    pred2 = Predicate('friends', 2)

    c1 = Constant('bob')
    v1 = LogicalVariable('X')
    atom1 = Atom(pred1, [c1])
    atom2 = Atom(pred2, [c1, v1])

    literal1 = Literal(atom1, True)
    literal2 = Literal(atom2, False)

    uclause = UnconstrainedClause([literal1, literal2])

    # v1 = LogicalVariable('X')
    # c1 = Constant('bob')
    c2 = Constant('a')
    d = SetOfConstants([c1, c2])

    eq_constraint = EqualityConstraint(v1, c2)
    ineq_constraint = InequalityConstraint(v1, c2)
    in_constraint = InclusionConstraint(c1, d)
    notin_constraint = NotInclusionConstraint(c1, d)

    cs = ConstraintSet([eq_constraint, ineq_constraint, in_constraint, notin_constraint])

    v2 = LogicalVariable('Y')
    clause = ConstrainedClause(uclause, [v1, v2], cs)
    print(clause)

    







