"""
Classes for clauses in FOL-DC.
This includes constrained AND unconstrained clauses.
"""

from logicalterms import *
from literals import *
from constraintsets import *
from abc import ABC

class Clause(ABC):
    """Abstract base class for constrained and unconstrained clauses"""

class UnconstrainedClause(Clause):
    """An FOL unconstrained clause.
    This consists of a set of FOL literals, which form a disjunction
    """

    def __init__(self, literals: list['Literal']) -> None:
        self.literals = literals

    def __str__(self) -> str:
        literal_strs = [str(literal) for literal in self.literals]
        return f"({' OR '.join(literal_strs)})"


class ConstrainedClause(Clause):
    """An FOL-DC constrained clause.
    This consists of an unconstrained clause with a set of bound variables and a constraint set.

    NOTE: For now we only work with bound *logical* variables."""

    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: list['LogicalVariable'],
            cs: 'ConstraintSet'):
        self.unconstrained_clause = unconstrained_clause
        self.bound_vars = bound_vars
        self.cs = cs

    def __str__(self) -> str:
        bound_vars_strs = [str(var) for var in self.bound_vars]
        return f"FORALL {{{', '.join(bound_vars_strs)}}}, {self.cs} : {self.unconstrained_clause}"



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

    







