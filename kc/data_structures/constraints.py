"""
Classes for domain constraints and constraint sets.

NOTE: at the moment only the types of constraints used in Chapter 4 of the PhD
are considered.
"""

from kc.data_structures.logicalterms import *
from kc.data_structures.domainterms import *
from abc import ABC

from typing import List

class ConstraintSet:
    """A FOL-DC constraint set.
    This consists of a set of constraints, which form a conjunction.
    """

    def __init__(self, constraints: List['Constraint']) -> None:
        self.constraints = constraints

    def __str__(self) -> str:
        constraint_strs = [f'({str(constraint)})' for constraint in self.constraints]
        return f"({' AND '.join(constraint_strs)})"

    def __repr__(self) -> str:
        return self.__str__()

class Constraint(ABC):
    """
    Abstract base class for constraints.
    This covers equality constraints, inclusion constraints, and their negations (plus more if needed).
    """

class LogicalConstraint(Constraint):
    """Abstract base class for constraints that only involve logical terms.
    This covers equality constraints and ineqality constraints."""

class SetConstraint(Constraint):
    """Abstract base class for constraints that involve domain terms (i.e. sets of constants and variables representing sets)
    This covers inclusion constraints and non-inclusion constraints."""

class EqualityConstraint(LogicalConstraint):
    """
    A FOL-DC equality constraint (between logical terms).
    This consists of two logical terms
    """

    def __init__(self, left_term: 'LogicalTerm', right_term: 'LogicalTerm') -> None:
        self.left_term = left_term
        self.right_term = right_term

    def __str__(self) -> str:
        return f'{self.left_term} == {self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

class InequalityConstraint(LogicalConstraint):
    """
    A FOL-DC inequality constraint (between logical terms).
    This consists of two logical terms
    """

    def __init__(self, left_term: 'LogicalTerm', right_term: 'LogicalTerm') -> None:
        self.left_term = left_term
        self.right_term = right_term

    def __str__(self) -> str:
        return f'{self.left_term} != {self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

class InclusionConstraint(SetConstraint):
    """
    A FOL-DC inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalTerm', domain_term: 'SetOfConstants') -> None:
        """NOTE: For now, we only allow 'SetOfConstants' rather than general 'DomainTerm'
        for the domain term"""
        self.logical_term = logical_term
        self.domain_term = domain_term

    def __str__(self) -> str:
        return f'{self.logical_term} IN {self.domain_term}'

    def __repr__(self) -> str:
        return self.__str__()

class NotInclusionConstraint(SetConstraint):
    """
    A FOL-DC negated inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalTerm', domain_term: 'SetOfConstants') -> None:
        """NOTE: For now, we only allow 'SetOfConstants' rather than general 'DomainTerm'
        for the domain term"""
        self.logical_term = logical_term
        self.domain_term = domain_term

    def __str__(self) -> str:
        return f'{self.logical_term} NOT IN {self.domain_term}'

    def __repr__(self) -> str:
        return self.__str__()

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




