"""
Classes for clauses in FOL-DC.
This includes constrained AND unconstrained clauses.
"""

from logicalterms import *
from literals import *
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

if __name__ == '__main__':
    pred1 = Predicate('smokes', 1)
    pred2 = Predicate('friends', 2)

    c1 = Constant('bob')
    v1 = LogicalVariable('X')
    atom1 = Atom(pred1, [c1])
    atom2 = Atom(pred2, [c1, v1])

    literal1 = Literal(atom1, True)
    literal2 = Literal(atom2, False)

    clause = UnconstrainedClause([literal1, literal2])
    print(clause)



