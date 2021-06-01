"""
Classes for domain constraints and constraint sets.

NOTE: at the moment only the types of constraints used in Chapter 4 of the PhD
are considered.
"""

from logicalterms import *
from domainterms import *
from abc import ABC

class ConstraintSet:
    """A FOL-DC constraint set.
    This consists of a set of constraints, which form a conjunction.
    """

    def __init__(self, constraints: list['Constraint']) -> None:
        self.constraints = constraints

    def __str__(self) -> str:
        constraint_strs = [f'({str(constraint)})' for constraint in self.constraints]
        return f"({' AND '.join(constraint_strs)})"


class Constraint(ABC):
    """
    Abstract base class for constraints.
    This covers equality constraints, inclusion constraints, and their negations (plus more if needed).
    """

class EqualityConstraint(Constraint):
    """
    A FOL-DC equality constraint (between logical terms).
    This consists of two logical terms
    """

    def __init__(self, left_term: 'LogicalTerm', right_term: 'LogicalTerm') -> None:
        self.left_term = left_term
        self.right_term = right_term

    def __str__(self) -> str:
        return f'{self.left_term} == {self.right_term}'


class InequalityConstraint(Constraint):
    """
    A FOL-DC inequality constraint (between logical terms).
    This consists of two logical terms
    """

    def __init__(self, left_term: 'LogicalTerm', right_term: 'LogicalTerm') -> None:
        self.left_term = left_term
        self.right_term = right_term

    def __str__(self) -> str:
        return f'{self.left_term} != {self.right_term}'


class InclusionConstraint(Constraint):
    """
    A FOL-DC inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalTerm', domain_term: 'DomainTerm') -> None:
        self.logical_term = logical_term
        self.domain_term = domain_term

    def __str__(self) -> str:
        return f'{self.logical_term} IN {self.domain_term}'


class NotInclusionConstraint(Constraint):
    """
    A FOL-DC negated inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalTerm', domain_term: 'DomainTerm') -> None:
        self.logical_term = logical_term
        self.domain_term = domain_term

    def __str__(self) -> str:
        return f'{self.logical_term} NOT IN {self.domain_term}'


if __name__ == '__main__':
    v1 = LogicalVariable('X')
    c1 = Constant('bob')
    c2 = Constant('a')

    eq_constraint = EqualityConstraint(v1, c2)
    print(eq_constraint)
    ineq_constraint = InequalityConstraint(v1, c2)
    print(ineq_constraint)

    d = SetOfConstants([c1, c2])
    in_constraint = InclusionConstraint(c1, d)
    print(in_constraint)

    notin_constraint = NotInclusionConstraint(c1, d)
    print(notin_constraint)

    cs = ConstraintSet([eq_constraint, ineq_constraint, in_constraint, notin_constraint])
    print(cs)




