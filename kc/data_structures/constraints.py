"""
Classes for domain constraints and constraint sets.

NOTE: at the moment only the types of constraints used in Chapter 4 of the PhD
are considered.
"""

from kc.data_structures import EquivalenceClass, EquivalenceClasses, LogicalVariable, SetOfConstants, Constant, VariableEquivalenceClass

from abc import ABC, abstractmethod

from typing import List, Iterable, Any, FrozenSet, Tuple, Union, Optional, Sequence

# defining type alias to simplify type hinting
ECList = List['EquivalenceClass']
VECSeq = Sequence['VariableEquivalenceClass']

class ConstraintSet:
    """A FOL-DC constraint set.
    This consists of a set of constraints, which form a conjunction.
    """

    def __init__(self, constraints: Iterable['Constraint']) -> None:
        self._constraints = frozenset(constraints)

    @property
    def constraints(self) -> FrozenSet['Constraint']:
        return self._constraints

    def join(self, other: 'ConstraintSet') -> 'ConstraintSet':
        """Create a ConstraintSet by joining together the constraints of two constraint sets"""
        new_constraints = self.constraints.union(other.constraints)
        return ConstraintSet(new_constraints)

    # def apply_substitution(self, substitution: 'Substitution') -> 'ConstraintSet':
    #     """Create a ConstraintSet by applying a Substitution to the constraints"""
    #     new_constraints: List['Constraint'] = []
    #     for constraint in self.constraints:
    #         new_constraint = constraint.apply_substitution(substitution)
    #         new_constraints.append(new_constraint)
    #     return ConstraintSet(new_constraints)

    def is_satisfiable(self) -> bool:
        """Is this constraint set satisfiable? i.e. are there any substitutions to
        its variables that do not contain contradictions?
        NOTE TODO: There are three major flaws with this function still:
            1) It cannot handle domain variables
            2) If there are more mutually unequal variables in a domain than there are elements
            of that domain it is not satisfiable but this function cannot determine that yet
            3) Inequality constraints between free variables aren't checked if they don't
            apeear in any equality constraints"""
        # First, we construct variable equivalence classes from the equality constraints 
        eq_classes = self.get_var_eq_classes()
        # Next, we make sure that these are consistent with the inequality constraints
        if not consistent_with_inequality_constraints(eq_classes, self):
            return False
        # Then we construct the possible domains for each equivalence class and
        # check that they are non-empty
        if not consistent_with_set_constraints(eq_classes, self):
            return False
        # finally, we assume that it is satisfiable
        return True

    def get_var_eq_classes(self) -> 'EquivalenceClasses[VariableEquivalenceClass]':
        """Return the variable equivalence classes given by the equality constraints
        of this constraint set
        NOTE: This uses the assumption that equality constraints only contain logical variables"""
        initial_eq_class_pairs: List['VariableEquivalenceClass'] = [] # collect the equivalence classes we will use later
        for constraint in self:
            if isinstance(constraint, EqualityConstraint):
                initial_eq_class_pairs.append(VariableEquivalenceClass(constraint.terms))

        initial_eq_classes: EquivalenceClasses['VariableEquivalenceClass'] = EquivalenceClasses(initial_eq_class_pairs)
        final_eq_classes = initial_eq_classes.make_self_consistent()
        return final_eq_classes

    def __iter__(self):
        return iter(self.constraints)

    def __eq__(self, other: Any) -> bool:
        """Two constraint sets are equal if their constraints are equal."""
        if not isinstance(other, ConstraintSet):
            return False
        same_constraints = (self.constraints == other.constraints)
        return same_constraints

    def __hash__(self) -> int:
        return hash(self.constraints)

    def __str__(self) -> str:
        constraint_strs = [f'({str(constraint)})' for constraint in self.constraints]
        logical_and_string = ' \u2227 '
        return f"({logical_and_string.join(constraint_strs)})"

    def __repr__(self) -> str:
        return self.__str__()

class Constraint(ABC):
    """
    Abstract base class for constraints.
    This covers equality constraints, inclusion constraints, and their negations (plus more if needed).
    """
    
    terms: Tuple['LogicalVariable', Union['LogicalVariable', 'SetOfConstants']] # every constraint needs some terms

    @abstractmethod
    def __hash__(self) -> int:
        """Constraints need to be hashable so we can put them in sets"""
        pass

    # @abstractmethod 
    # def apply_substitution(self, substitution: 'Substitution') -> 'Constraint':
    #     """We need to be able to apply substitutions to constraints"""
    #     pass

    @abstractmethod
    def contains_contradiction(self) -> bool:
        """Does this constraint contain a contradiction?
        TODO: decide if this should be a property"""
        pass


class LogicalConstraint(Constraint):
    """Abstract base class for constraints that only involve logical terms.
    This covers equality constraints and ineqality constraints.
    NOTE: I am currently experimenting with only having logical variables in logical constraints,
    no constants."""

    @property
    @abstractmethod
    def left_term(self) -> 'LogicalVariable':
        pass

    @property
    @abstractmethod
    def right_term(self) -> 'LogicalVariable':
        pass

    def __hash__(self) -> int:
        return hash((self.left_term, self.right_term))


class SetConstraint(Constraint):
    """Abstract base class for constraints that involve domain terms (i.e. sets of constants and variables representing sets)
    This covers inclusion constraints and non-inclusion constraints.

    NOTE TODO: For now, this doesn't allow domain variables"""

    @property
    @abstractmethod
    def logical_term(self) -> 'LogicalVariable':
        pass

    @property
    @abstractmethod
    def domain_term(self) -> 'SetOfConstants':
        pass

    def __hash__(self) -> int:
        """TODO: see if this makes sense as a hash function"""
        return hash((self.logical_term, self.domain_term))


class EqualityConstraint(LogicalConstraint):
    """
    A FOL-DC equality constraint (between logical terms).
    This consists of two logical terms
    """

    def __init__(self, left_term: 'LogicalVariable', right_term: 'LogicalVariable') -> None:
        self.terms: Tuple['LogicalVariable', 'LogicalVariable'] = (left_term, right_term)

    @property
    def left_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def right_term(self) -> 'LogicalVariable':
        return self.terms[1]

    # def apply_substitution(self, substitution: 'Substitution') -> 'EqualityConstraint':
    #     """Apply a substitution to the constraint, returning a new constraint.
    #     TODO: Figure out how to abstract this up the the base class"""
    #     return EqualityConstraint(substitution[self.left_term], substitution[self.right_term])

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For EqualityConstraint, this means checking if the two sides are different constants"""
        both_terms_constants = isinstance(self.left_term, Constant) and isinstance(self.right_term, Constant)
        if both_terms_constants and self.left_term != self.right_term:
            return True
        else:
            return False

    def __eq__(self, other: Any) -> bool:
        """Two equality constraints are equal if they mention the same terms (note the order doesn't matter)"""
        if not isinstance(other, EqualityConstraint):
            return False
        same_way = (self.left_term == other.left_term and self.right_term == other.right_term)
        flipped = (self.left_term == other.right_term and self.right_term == other.left_term)
        return same_way or flipped

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between EqualityConstraints and InequalityConstraints"""
        return super().__hash__()

    def __invert__(self) -> 'InequalityConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !="""
        return InequalityConstraint(self.left_term, self.right_term)

    def __str__(self) -> str:
        return f'{self.left_term} = {self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

class InequalityConstraint(LogicalConstraint):
    """
    A FOL-DC inequality constraint (between logical terms).
    This consists of two logical terms
    """

    def __init__(self, left_term: 'LogicalVariable', right_term: 'LogicalVariable') -> None:
        self.terms: Tuple['LogicalVariable', 'LogicalVariable'] = (left_term, right_term)

    @property
    def left_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def right_term(self) -> 'LogicalVariable':
        return self.terms[1]

    # def apply_substitution(self, substitution: 'Substitution') -> 'InequalityConstraint':
    #     """Apply a substitution to the constraint, returning a new constraint.
    #     TODO: Figure out how to abstract this up the the base class"""
    #     return InequalityConstraint(substitution[self.left_term], substitution[self.right_term])

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For InequalityConstraint, this means checking if the two sides are the same term"""
        return self.left_term == self.right_term

    def __eq__(self, other: Any) -> bool:
        """Two inequality constraints are equal if they mention the same terms (note the order doesn't matter)"""
        if not isinstance(other, InequalityConstraint):
            return False
        same_way = (self.left_term == other.left_term and self.right_term == other.right_term)
        flipped = (self.left_term == other.right_term and self.right_term == other.left_term)
        return same_way or flipped

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between EqualityConstraints and InequalityConstraints"""
        return super().__hash__()

    def __invert__(self) -> 'EqualityConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !="""
        return EqualityConstraint(self.left_term, self.right_term)

    def __str__(self) -> str:
        not_equal_string = ' \u2260 '
        return f'{self.left_term}{not_equal_string}{self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

class InclusionConstraint(SetConstraint):
    """
    A FOL-DC inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalVariable', domain_term: 'SetOfConstants') -> None:
        """NOTE TODO: For now, we only allow 'SetOfConstants' rather than general 'DomainTerm'
        for the domain term"""
        self.terms: Tuple['LogicalVariable', 'SetOfConstants'] = (logical_term, domain_term)

    @property
    def logical_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def domain_term(self) -> 'SetOfConstants':
        return self.terms[1]

    # def apply_substitution(self, substitution: 'Substitution') -> 'InclusionConstraint':
    #     """Apply a substitution to the constraint, returning a new constraint.
    #     TODO: Figure out how to abstract this up the the base class.
    #     TODO: Include domain substitutions"""
        # return InclusionConstraint(substitution[self.logical_term], self.domain_term)

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For InclusionConstraint, this means checking if the logical term is a constant and the domain term a set of constants without the constant"""
        if isinstance(self.logical_term, Constant) and isinstance(self.domain_term, SetOfConstants):
            return not self.logical_term in self.domain_term
        else:
            return False

    def __eq__(self, other: Any) -> bool:
        """Two inclusion constraints are equal if they mention the same logical and domain terms """
        if not isinstance(other, InclusionConstraint):
            return False
        return self.logical_term == other.logical_term and self.domain_term == other.domain_term

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between InclusionConstraints and NotInclusionConstraints"""
        return super().__hash__()

    def __invert__(self) -> 'NotInclusionConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !="""
        return NotInclusionConstraint(self.logical_term, self.domain_term)

    def __str__(self) -> str:
        element_of_string = ' \u2208 '
        return f'{self.logical_term}{element_of_string}{self.domain_term}'

    def __repr__(self) -> str:
        return self.__str__()

class NotInclusionConstraint(SetConstraint):
    """
    A FOL-DC negated inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalVariable', domain_term: 'SetOfConstants') -> None:
        """NOTE: For now, we only allow 'SetOfConstants' rather than general 'DomainTerm'
        for the domain term"""
        self.terms: Tuple['LogicalVariable', 'SetOfConstants'] = (logical_term, domain_term)

    @property
    def logical_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def domain_term(self) -> 'SetOfConstants':
        return self.terms[1]

    # def apply_substitution(self, substitution: 'Substitution') -> 'NotInclusionConstraint':
    #     """Apply a substitution to the constraint, returning a new constraint.
    #     TODO: Figure out how to abstract this up the the base class.
    #     TODO: Include domain substitutions"""
    #     return NotInclusionConstraint(substitution[self.logical_term], self.domain_term)

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For NotInclusionConstraint, this means checking if the logical term is a constant and the domain term a set of constants containing the constant"""
        if isinstance(self.logical_term, Constant) and isinstance(self.domain_term, SetOfConstants):
            return self.logical_term in self.domain_term
        else:
            return False

    def __eq__(self, other: Any) -> bool:
        """Two not-inclusion constraints are equal if they mention the same logical and domain terms """
        if not isinstance(other, NotInclusionConstraint):
            return False
        return self.logical_term == other.logical_term and self.domain_term == other.domain_term

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between InclusionConstraints and NotInclusionConstraints"""
        return super().__hash__()

    def __invert__(self) -> 'InclusionConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !="""
        return InclusionConstraint(self.logical_term, self.domain_term)

    def __str__(self) -> str:
        not_element_of_string = ' \u2209 '
        return f'{self.logical_term}{not_element_of_string}{self.domain_term}'

    def __repr__(self) -> str:
        return self.__str__()
