"""
Class for FO-CNF formulas.
"""

from kc.data_structures.clauses import *

from typing import List

class FO_CNF:
    """
    A FOL-DC CNF.
    This consists of a set of constrained clauses, which form a conjunction.
    """

    def __init__(self, clauses: List['ConstrainedClause']) -> None:
        self.clauses = clauses

    def __str__(self) -> str:
        clause_strs = [f'({str(clause)})' for clause in self.clauses]
        return '\nAND\n'.join(clause_strs)


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

    cs1 = ConstraintSet([eq_constraint, in_constraint])
    cs2 = ConstraintSet([ineq_constraint, notin_constraint])

    v2 = LogicalVariable('Y')
    v3 = LogicalVariable('Z')
    clause1 = ConstrainedClause(uclause, [v1, v2], cs1)
    clause2 = ConstrainedClause(uclause, [v2, v3], cs2)

    cnf = FO_CNF([clause1, clause2])
    print(cnf)
    







