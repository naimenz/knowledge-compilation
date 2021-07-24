"""
Classes for domain constraints and constraint sets AND NOW Substitutions, because they depend on each other. 
AND NOW ALSO EQUIVALENCE CLASSES
NOTE: at the moment only the types of constraints used in Chapter 4 of the PhD
are considered.
"""

from kc.data_structures import LogicalVariable, SetOfConstants, Constant, DomainTerm, DomainVariable, ProperDomain, EmptyDomain
from kc.util import powerset, get_element_of_set

from abc import ABC, abstractmethod
from functools import reduce

from typing import List, Iterable, Any, FrozenSet, Tuple, Union, Optional, Sequence, Set
from typing import TypeVar, cast, Generic
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import CNF, LogicalVariable, LogicalTerm, ConstrainedAtom

VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']
TEC = TypeVar('TEC', bound='EquivalenceClass') # this is a type var so I can work with eq classes or var eq classes

class ConstraintSet:
    """A FOL-DC constraint set.
    This consists of a set of constraints, which form a conjunction.
    """

    def __init__(self, constraints: Iterable['Constraint']) -> None:
        # make it easier to access different types of constraints 
        equality_constraints, inequality_constraints = set(), set()
        inclusion_constraints, notinclusion_constraints = set(), set()

        for constraint in constraints:
            if isinstance(constraint, EqualityConstraint):
                equality_constraints.add(constraint)
            elif isinstance(constraint, InequalityConstraint):
                inequality_constraints.add(constraint)
            elif isinstance(constraint, InclusionConstraint):
                inclusion_constraints.add(constraint)
            elif isinstance(constraint, NotInclusionConstraint):
                notinclusion_constraints.add(constraint)

        redundant_inclusion_constraints = self._get_redundant_inclusion_constraints(inclusion_constraints)

        self._inequality_constraints = frozenset(inequality_constraints)
        self._equality_constraints = frozenset(equality_constraints)
        self._inclusion_constraints = frozenset(inclusion_constraints) - redundant_inclusion_constraints
        self._notinclusion_constraints = frozenset(notinclusion_constraints)

        self._constraints = frozenset(constraints) - redundant_inclusion_constraints
    
    def _get_redundant_inclusion_constraints(self,
                                             inclusion_constraints: Set['InclusionConstraint']
                                             ) -> Set['InclusionConstraint']:
        """Get the inclusion constraints that are already implied by some stricter constraint"""
        redundant_inclusion_constraints = set()
        for ic in inclusion_constraints:
            if isinstance(ic.logical_term, LogicalVariable) and isinstance(ic.domain_term, ProperDomain):
                variable = ic.logical_term
                domain = ic.domain_term
                for other_ic in inclusion_constraints:
                    if other_ic.logical_term == variable \
                            and isinstance(other_ic.domain_term, ProperDomain) \
                            and domain.is_strict_superset_of(other_ic.domain_term):
                        redundant_inclusion_constraints.add(ic)
                        break
        return redundant_inclusion_constraints
    

    @property
    def constraints(self) -> FrozenSet['Constraint']:
        return self._constraints

    @property
    def equality_constraints(self) -> FrozenSet['EqualityConstraint']:
        return self._equality_constraints

    @property
    def inequality_constraints(self) -> FrozenSet['InequalityConstraint']:
        return self._inequality_constraints

    @property
    def inclusion_constraints(self) -> FrozenSet['InclusionConstraint']:
        return self._inclusion_constraints

    @property
    def notinclusion_constraints(self) -> FrozenSet['NotInclusionConstraint']:
        return self._notinclusion_constraints

    @property
    def logical_constraints(self) -> FrozenSet['LogicalConstraint']:
        eq = cast(FrozenSet['LogicalConstraint'], self.equality_constraints) # hack for type checking
        ineq = cast(FrozenSet['LogicalConstraint'], self.inequality_constraints) # hack for type checking
        return eq.union(ineq)

    @property
    def set_constraints(self) -> FrozenSet['SetConstraint']:
        inc = cast(FrozenSet['SetConstraint'], self.inclusion_constraints) # hack for type checking
        notinc = cast(FrozenSet['SetConstraint'], self.notinclusion_constraints) # hack for type checking
        return inc.union(notinc)


    def unequal_constants_for_variable(self, variable: 'LogicalVariable') -> Set['Constant']:
        """Return all the constants that 'variable' is not equal to"""
        unequal_constants: Set['Constant'] = set()
        for constraint in self.notinclusion_constraints:
            if constraint.logical_term == variable and isinstance(constraint.domain_term, SetOfConstants):
                unequal_constants.update(constraint.domain_term.constants)
        return unequal_constants


    def join(self, other: 'ConstraintSet') -> 'ConstraintSet':
        """Create a ConstraintSet by joining together the constraints of two constraint sets"""
        new_constraints = self.constraints.union(other.constraints)
        return ConstraintSet(new_constraints)

    def is_non_empty(self) -> bool:
        return len(self.constraints) > 0

    def substitute(self, substitution: 'Substitution') -> Optional['ConstraintSet']:
        """Create a ConstraintSet by applying a Substitution to the constraints
        NOTE: If after substituting the cs is unsatisfiable, return None
        NOTE: this may be counter-intuitive because applying a substitution
        to an EqualityConstraint may produce an InclusionConstraint"""
        new_constraints: List['Constraint'] = []
        for constraint in self.constraints:
            new_constraint = constraint.substitute(substitution)
            new_constraints.append(new_constraint)

        if FalseConstraint() in new_constraints:
            print(f'Substitution {substitution} made cs unsatisfiable!')
            return None
            # raise ValueError(f'Substitution {substitution} made cs unsatisfiable!')

        # remove trivial constraints (always true)
        filtered_constraints = [c for c in new_constraints if c != EmptyConstraint()]

        return ConstraintSet(filtered_constraints)

    def drop_constraints_involving_only_specific_variables(self, variables: Iterable['LogicalVariable']) -> 'ConstraintSet':
        """Return a constraint set where all constraints that only involve variables from 'variables' are removed"""
        relevant_logical_constraints = [lc for lc in self.logical_constraints if (not lc.left_term in variables and not lc.right_term in variables)]
        relevant_set_constraints = [sc for sc in self.set_constraints if not sc.logical_term in variables]
        return ConstraintSet([*relevant_logical_constraints, *relevant_set_constraints])

    def project(self, c_atom: 'ConstrainedAtom') -> 'ConstraintSet':
        """Return a constraint set with only the constraints involving a specific c_atom"""
        relevant_variables = c_atom.literal_variables  # TODO: replace with non-property version
        relevant_logical_constraints = [lc for lc in self.logical_constraints if (lc.left_term in relevant_variables and lc.right_term in relevant_variables)]
        relevant_set_constraints = [sc for sc in self.set_constraints if sc.logical_term in relevant_variables]
        return ConstraintSet([*relevant_logical_constraints, *relevant_set_constraints])

    def get_domain_for_variable(self, variable: 'LogicalVariable') -> 'ProperDomain':
        """Get the domain for a specific variable in this constraint set"""
        domains = []
        for constraint in self.inclusion_constraints:
            if constraint.logical_term == variable and isinstance(constraint.domain_term, ProperDomain):
                domains.append(constraint.domain_term)
        return ProperDomain.intersect_all(*domains)
        raise ValueError(f"{variable} has no ProperDomain")

    def get_allowed_constants_for(self, variable: 'LogicalVariable') -> FrozenSet['Constant']:
        """Get the constants that this variable could be equal to in this constraint set
        NOTE: This struggles with DomainVariables and partially with free variables"""
        included_domains = []
        excluded_domains = []
        for constraint in self:
            if isinstance(constraint, InclusionConstraint) and constraint.logical_term == variable:
                included_domains.append(constraint.domain_term)
            elif isinstance(constraint, NotInclusionConstraint) and constraint.logical_term == variable:
                excluded_domains.append(constraint.domain_term)
        included_domain = DomainTerm.intersect_constants(*included_domains)
        excluded_domain = DomainTerm.union_constants(*excluded_domains)
        variable_domain = included_domain.difference(excluded_domain)
        return variable_domain


    def is_satisfiable(self) -> bool:
        """Is this constraint set satisfiable? i.e. are there any substitutions to
        its variables that do not contain contradictions?
        NOTE TODO: Redo this function altogether? I'm not sure how to approach this"""
        if any(isinstance(constraint, FalseConstraint) for constraint in self.constraints):
            # DEBUG
            fcs = [c for c in self.constraints if isinstance(c, FalseConstraint)]
            print([fc.debug_message for fc in fcs])
            return False

        # a bit of a hack but can catch some cases
        for c in self.constraints:
            if ~c in self.constraints:
                return False

        eq_classes = self.get_var_eq_classes()
        if not eq_classes.consistent_with_variable_inequality_constraints(self):
            return False

        if not eq_classes.consistent_with_constant_inequality_constraints(self):
            return False

        # # Then we construct the possible domains for each equivalence class and
        # # check that they are non-empty
        if not eq_classes.consistent_with_set_constraints(self):
            return False

        if self._too_many_mutually_unequal():
            return False
        return True

    def _too_many_mutually_unequal(self) -> bool:
        """ This is a check that we don't have a situation like:
        X in {a}, Y in {a}, X != Y
        i.e. that the domain for unequal variables isn't too small to support inequalities
        between all of the variables with that domain
        NOTE: we only check variables with an InclusionConstraint"""
        variables_with_domains = set(c.logical_term for c in self.set_constraints if isinstance(c, InclusionConstraint))
        inequality_constraints = set(c for c in self.logical_constraints if isinstance(c, InequalityConstraint))
        mutual_inequalities = self._find_mutual_inequalities(variables_with_domains, inequality_constraints)

        # we check for too many variables in too small a domain by taking the union of the domains
        for mutual_inequality in mutual_inequalities:
            domains_generator = (self.get_allowed_constants_for(variable) for variable in mutual_inequality)
            combined_domain = frozenset.union(*domains_generator)
            if len(combined_domain) < len(mutual_inequality):
                return True
        return False

    def _find_mutual_inequalities(self, variables: Set['LogicalVariable'],
                                  inequalities: Set['InequalityConstraint']
                                  ) -> Set[FrozenSet['LogicalVariable']]:
        """Given a set of inequalities and variables that we care about,
         find every set of mutual inequalities (e.g. X != Y, Y != Z, X != Z)"""
        mutual_inequalities = set()
        for subset in powerset(variables):
            if len(subset) < 2:
                continue
            if self._subset_are_all_mutually_unequal(subset, inequalities):
                mutual_inequalities.add(frozenset(subset))
        return mutual_inequalities


    def _subset_are_all_mutually_unequal(self, subset: Iterable['LogicalVariable'],
                                         inequalities: Set['InequalityConstraint']) -> bool:
        """Check if a given subset of variables is mutually unequal"""
        # surprisingly, two for-loops are faster than using itertools.product
        for variable in subset:
            for other_variable in subset:
                if variable != other_variable and not ( InequalityConstraint(variable, other_variable) in inequalities ):
                    return False
        return True


    def get_var_eq_classes(self) -> 'EquivalenceClasses[VariableEquivalenceClass]':
        """Return the variable equivalence classes given by the equality constraints
        of this constraint set (including singletons).
        NOTE: This uses the assumption that equality constraints only contain logical variables"""
        initial_eq_class_pairs: List['VariableEquivalenceClass'] = [] # collect the equivalence classes we will use later
        # first add a class for each individual variable
        for var in self.get_logical_variables():
            initial_eq_class_pairs.append(VariableEquivalenceClass((var, var)))

        for constraint in self:
            if isinstance(constraint, EqualityConstraint):
                initial_eq_class_pairs.append(VariableEquivalenceClass(constraint.terms))

        initial_eq_classes: EquivalenceClasses['VariableEquivalenceClass'] = EquivalenceClasses(initial_eq_class_pairs)
        final_eq_classes = initial_eq_classes.make_self_consistent()
        return final_eq_classes


    def get_logical_variables(self) -> Set['LogicalVariable']:
        """Extract just the variables from each constraint in the constraint set"""
        logical_variables: Set['LogicalVariable'] = set()
        for constraint in self:
            for term in constraint.terms:
                if isinstance(term, LogicalVariable):
                    logical_variables.add(term)
        return logical_variables

    def __iter__(self):
        return iter(self.constraints)

    def __eq__(self, other: Any) -> bool:
        """Two constraint sets are equal if their constraints are equal."""
        return isinstance(other, ConstraintSet) and self.constraints == other.constraints

    def __hash__(self) -> int:
        return hash(self.constraints)

    def __str__(self) -> str:
        constraint_strs = [f'({str(constraint)})' for constraint in sorted(self.constraints)]
        logical_and_string = ' \u2227 '
        return f"({logical_and_string.join(constraint_strs)})"

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """Order of constraint set must depend on the orders of its constraints"""
        if not isinstance(other, ConstraintSet):
            raise NotImplementedError(f'Cannot compare ConstraintSet with {type(other)}')
        return sorted(self.constraints) < sorted(other.constraints)

class Constraint(ABC):
    """
    Abstract base class for constraints.
    This covers equality constraints, inclusion constraints, and their negations (plus more if needed).
    """
    
    terms: Tuple[Union['LogicalTerm', 'DomainTerm'], Union['LogicalTerm', 'DomainTerm']] 

    @abstractmethod
    def __hash__(self) -> int:
        """Constraints need to be hashable so we can put them in sets"""
        pass

    @abstractmethod 
    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """We need to be able to apply substitutions to constraints"""
        pass

    @abstractmethod
    def contains_contradiction(self) -> bool:
        """Does this constraint contain a contradiction?
        TODO: decide if this should be a property"""
        pass

    @abstractmethod
    def __invert__(self) -> 'Constraint':
        """Unfortunately we cannot be more specific about the return type because 
        each constraint returns its opposite"""
        pass

    @abstractmethod
    def __lt__(self, other: Any) -> bool:
        """All constraints need to be comparable"""
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
        """We hash it as a set with the class to make sure
        flipped constraints have the same hash, and Equality/InequalityConstraint
        have different hashes"""
        return hash((self.__class__, frozenset((self.left_term, self.right_term))))


class SetConstraint(Constraint):
    """Abstract base class for constraints that involve domain terms (i.e. sets of constants and variables representing sets)
    This covers inclusion constraints and non-inclusion constraints.

    NOTE TODO: For now, this doesn't work with domain variables"""

    @property
    @abstractmethod
    def logical_term(self) -> 'LogicalVariable':
        pass

    @property
    @abstractmethod
    def domain_term(self) -> 'DomainTerm':
        pass

    def __hash__(self) -> int:
        return hash((self.__class__, self.logical_term, self.domain_term))


class DomainConstraint(Constraint):
    """Abstract base class for constraints that involve a relationship between two domain terms.
    For now, this is just subset constraints"""

    @property
    @abstractmethod
    def left_term(self) -> 'DomainTerm':
        pass

    @property
    @abstractmethod
    def right_term(self) -> 'DomainTerm':
        pass

    def __hash__(self) -> int:
        """We hash it as a set with the class to make sure
        flipped constraints have the same hash, and different types of DomainConstraint 
        have different hashes"""
        return hash((self.__class__, frozenset((self.left_term, self.right_term))))

class SubsetConstraint(DomainConstraint):
    """A class for one domain term being a subset of another
    NOTE: For now I will enforce that the left term is a DomainVariable, but I don't know
    if this is strictly necessary"""

    def __init__(self, left_term: 'DomainVariable', right_term: 'DomainTerm') -> None:
        self.terms: Tuple['DomainVariable', 'DomainTerm'] = (left_term, right_term)

    @property
    def left_term(self) -> 'DomainVariable':
        return self.terms[0]

    @property
    def right_term(self) -> 'DomainTerm':
        return self.terms[1]

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        TODO: Include domain substitutions"""
        raise NotImplementedError("Hopefully won't substitute domains")

    def contains_contradiction(self) -> bool:
        raise NotImplementedError("Hopefully we won't be checking Subset constraints for contradictions")

    def __eq__(self, other: Any) -> bool:
        """Two inclusion constraints are equal if they mention the same logical and domain terms """
        return isinstance(other, SubsetConstraint) \
               and self.left_term == other.left_term \
               and self.right_term == other.right_term

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between InclusionConstraints and NotInclusionConstraints"""
        return super().__hash__()

    def __invert__(self) -> 'NotSubsetConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !="""
        return NotSubsetConstraint(self.left_term, self.right_term)

    def __str__(self) -> str:
        subset_of_string = ' \u2286 '
        return f'{self.left_term}{subset_of_string}{self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """Subset constraints appear first"""
        if not isinstance(other, Constraint):
            raise NotImplementedError(f'Cannot compare SubsetConstraint with {type(other)}')
        elif isinstance(other, LogicalConstraint) or isinstance(other, SetConstraint):
            return True
        elif isinstance(other, NotSubsetConstraint):
            return True
        else:
            return self.terms < other.terms

class NotSubsetConstraint(DomainConstraint):
    """A class for one domain term NOT being a subset of another
    NOTE: For now I will enforce that the left term is a DomainVariable, but I don't know
    if this is strictly necessary"""

    def __init__(self, left_term: 'DomainVariable', right_term: 'DomainTerm') -> None:
        self.terms: Tuple['DomainVariable', 'DomainTerm'] = (left_term, right_term)

    @property
    def left_term(self) -> 'DomainVariable':
        return self.terms[0]

    @property
    def right_term(self) -> 'DomainTerm':
        return self.terms[1]

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        TODO: Include domain substitutions"""
        raise NotImplementedError("Hopefully won't substitute domains")

    def contains_contradiction(self) -> bool:
        raise NotImplementedError("Hopefully we won't be checking Subset constraints for contradictions")

    def __eq__(self, other: Any) -> bool:
        """Two inclusion constraints are equal if they mention the same logical and domain terms """
        return isinstance(other, NotSubsetConstraint) \
               and self.left_term == other.left_term \
               and self.right_term == other.right_term

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.
        """
        return super().__hash__()

    def __invert__(self) -> 'SubsetConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !="""
        return SubsetConstraint(self.left_term, self.right_term)

    def __str__(self) -> str:
        not_subset_of_string = ' \u2288 '
        return f'{self.left_term}{not_subset_of_string}{self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """Subset constraints appear first"""
        if not isinstance(other, Constraint):
            raise NotImplementedError(f'Cannot compare NotSubsetConstraint with {type(other)}')
        elif isinstance(other, LogicalConstraint) or isinstance(other, SetConstraint):
            return True
        elif isinstance(other, SubsetConstraint):
            return False
        else:
            return self.terms < other.terms

class EmptyConstraint(Constraint):
    """This is a special class for a constraint that is trivially satisfied (i.e. always true)"""
    def __init__(self, debug_message: str ='') -> None:
        """If desired, we can include a message explaining why this EmptyConstraint was made
        (for debugging)"""
        self.debug_message = debug_message

    def substitute(self, substitution: 'Substitution') -> 'EmptyConstraint':
        """Substitution does nothing because no terms"""
        return EmptyConstraint()
    
    def contains_contradiction(self) -> bool:
        """Never contains a contradiction because it's always true"""
        return False

    def __eq__(self, other: Any) -> bool:
        """All EmptyConstraints are equal, because the debug message is not important"""
        return isinstance(other, EmptyConstraint)

    def __hash__(self) -> int:
        """All EmptyConstraints are the same"""
        return hash('EmptyConstraint')

    def __str__(self) -> str:
        return f'EmptyConstraint({self.debug_message})'

    def __repr__(self) -> str:
        return self.__str__()

    def __invert__(self) -> 'FalseConstraint':
        return FalseConstraint('Inverted EmptyConstraint')

    def __lt__(self, other: Any) -> bool:
        """Shouldn't need < for this class"""
        raise NotImplementedError("Shouldn't need __lt__ for {self.__class__}")


class FalseConstraint(Constraint):
    """This is a special class for a constraint that is trivially UNsatisfied (i.e. always false)"""
    def __init__(self, debug_message: str ='') -> None:
        """If desired, we can include a message explaining why this FalseConstraint was made
        (for debugging)"""
        self.debug_message = debug_message

    def substitute(self, substitution: 'Substitution') -> 'FalseConstraint':
        """Substitution does nothing because no terms"""
        return FalseConstraint()
    
    def contains_contradiction(self) -> bool:
        """Always contains a contradiction because it's always false"""
        return True

    def __eq__(self, other: Any) -> bool:
        """Any two FalseConstraints are effectively equal.
        They could differ in debug message but this is not important"""
        return isinstance(other, FalseConstraint)

    def __hash__(self) -> int:
        """All EmptyConstraints are the same"""
        return hash('FalseConstraint')

    def __str__(self) -> str:
        return f'FalseConstraint({self.debug_message})'

    def __repr__(self) -> str:
        return self.__str__()

    def __invert__(self) -> 'EmptyConstraint':
        return EmptyConstraint('Inverted FalseConstraint')

    def __lt__(self, other: Any) -> bool:
        """Shouldn't need < for this class"""
        raise NotImplementedError("Shouldn't need __lt__ for {self.__class__}")

class EqualityConstraint(LogicalConstraint):
    """
    A FOL-DC equality constraint (between logical variables).
    This consists of two logical VARIABLES, not constants, which are handled by InclusionConstraint.
    """

    def __init__(self, left_term: 'LogicalVariable', right_term: 'LogicalVariable') -> None:
        self.terms: Tuple['LogicalVariable', 'LogicalVariable'] = (left_term, right_term)

    @property
    def left_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def right_term(self) -> 'LogicalVariable':
        return self.terms[1]

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        NOTE: This will return an InclusionConstraint if ONE of the new terms is a constant.
        TODO: Figure out how to abstract this up the the base class"""
        new_left_term = substitution[self.left_term]
        new_right_term = substitution[self.right_term]
        # avoid trivial constraints
        if new_left_term == new_right_term:
            return EmptyConstraint(f'{new_left_term} == {new_right_term}')

        if isinstance(new_left_term, Constant):
            if isinstance(new_right_term, Constant):
                if new_left_term != new_right_term:
                    return FalseConstraint(f'{new_left_term} != {new_right_term}')
                else:
                    # this branch won't be reached, just for type-checking
                    return EmptyConstraint(f'{new_left_term} == {new_right_term}')  
            else:
                right_var = cast('LogicalVariable', new_right_term) # hack for type checking
                return InclusionConstraint(right_var, SetOfConstants([new_left_term]))
        elif isinstance(new_right_term, Constant):
            left_var = cast('LogicalVariable', new_left_term) # hack for type checking
            return InclusionConstraint(left_var, SetOfConstants([new_right_term]))
        else:
            left_var = cast('LogicalVariable', new_left_term) # hack for type checking
            right_var = cast('LogicalVariable', new_right_term) # hack for type checking
            return EqualityConstraint(left_var, right_var)

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
        return isinstance(other, EqualityConstraint) \
               and ( (self.left_term == other.left_term and self.right_term == other.right_term) \
               or (self.left_term == other.right_term and self.right_term == other.left_term) )

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

    def __lt__(self, other: Any) -> bool:
        """Subset constraints appear first"""
        if isinstance(other, DomainConstraint) or isinstance(other, SetConstraint):
            return False
        elif isinstance(other, InequalityConstraint):
            return True
        elif isinstance(other, LessThanConstraint):
            return False
        elif isinstance(other, EqualityConstraint):
            return sorted(self.terms) < sorted(other.terms)
        else:
            raise NotImplementedError(f'Cannot compare EqualityConstraint with {type(other)}')

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

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        NOTE: This will return a NotInclusionConstraint if ONE of the new terms is a constant.
        """
        new_left_term = substitution[self.left_term]
        new_right_term = substitution[self.right_term]
        if isinstance(new_left_term, Constant):
            if isinstance(new_right_term, Constant):
                if new_left_term != new_right_term:
                    return EmptyConstraint(f'{new_left_term} != {new_right_term}')
                else:
                    return FalseConstraint(f'{new_left_term} == {new_right_term}')
            else:
                right_var = cast('LogicalVariable', new_right_term) # hack for type checking
                return NotInclusionConstraint(right_var, SetOfConstants([new_left_term]))
        elif isinstance(new_right_term, Constant):
            left_var = cast('LogicalVariable', new_left_term) # hack for type checking
            return NotInclusionConstraint(left_var, SetOfConstants([new_right_term]))
        else:
            left_var = cast('LogicalVariable', new_left_term) # hack for type checking
            right_var = cast('LogicalVariable', new_right_term) # hack for type checking
            return InequalityConstraint(left_var, right_var)

    def is_not_trivial(self, c_atom: 'ConstrainedAtom') -> bool:
        """Is this inequality constraint between variables 'trivial' given the domains
        of its variables? I.e. the constraint is trivial if the domains of the variables do not overlap,
        so they could never be equal anyway"""
        left_term_class = EquivalenceClass([self.left_term])
        right_term_class = EquivalenceClass([self.right_term])
        left_domain = left_term_class.get_shared_domain_from_cs(c_atom.cs)
        right_domain = right_term_class.get_shared_domain_from_cs(c_atom.cs)
        domain_terms_intersect = left_domain.intersect_with(right_domain).size() > 0
        return domain_terms_intersect

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For InequalityConstraint, this means checking if the two sides are the same term"""
        return self.left_term == self.right_term

    def __eq__(self, other: Any) -> bool:
        """Two inequality constraints are equal if they mention the same terms (note the order doesn't matter)"""
        return isinstance(other, InequalityConstraint) \
               and ( (self.left_term == other.left_term and self.right_term == other.right_term) \
               or (self.left_term == other.right_term and self.right_term == other.left_term) )

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

    def __lt__(self, other: Any) -> bool:
        """Subset constraints appear first"""
        if isinstance(other, DomainConstraint) or isinstance(other, SetConstraint):
            return False
        elif isinstance(other, EqualityConstraint):
            return False
        elif isinstance(other, LessThanConstraint):
            return False
        elif isinstance(other, InequalityConstraint):
            return sorted(self.terms) < sorted(other.terms)
        else:
            raise NotImplementedError(f'Cannot compare InequalityConstraint with {type(other)}')

class LessThanConstraint(LogicalConstraint):
    """This is a constraint on two variables (e.g. Y, Z) enforcing
    that Y < Z. This is interpreted as the natural order on the constants --
    since the set of constants is finite, we can always pick an order.
    NOTE: This constraint is ONLY used in the constraint sets of NNF nodes,
    so doesn't need any functionality"""

    def __init__(self, left_term: 'LogicalVariable', right_term: 'LogicalVariable') -> None:
        self.terms: Tuple['LogicalVariable', 'LogicalVariable'] = (left_term, right_term)

    @property
    def left_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def right_term(self) -> 'LogicalVariable':
        return self.terms[1]

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        """
        raise NotImplementedError('LessThanConstraint should not actually be used - it is for NNFNodes')

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For InequalityConstraint, this means checking if the two sides are the same term"""
        raise NotImplementedError('LessThanConstraint should not actually be used - it is for NNFNodes')


    def __eq__(self, other: Any) -> bool:
        """Two less-than constraints are equal if they mention the same terms (note the order doesn't matter)"""
        return isinstance(other, LessThanConstraint) \
               and ( (self.left_term == other.left_term and self.right_term == other.right_term) \
               or (self.left_term == other.right_term and self.right_term == other.left_term) )

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between EqualityConstraints and InequalityConstraints"""
        return super().__hash__()

    def __str__(self) -> str:
        less_than_string = ' < '
        return f'{self.left_term}{less_than_string}{self.right_term}'

    def __repr__(self) -> str:
        return self.__str__()

    def __invert__(self) -> 'Constraint':
        """Shouldn't need to invert LessThanConstraints because they are only
        really used in NNFnodes."""
        raise NotImplementedError('Cannot invert LessThanConstraint')

    def __lt__(self, other: Any) -> bool:
        """Subset constraints appear first"""
        if isinstance(other, DomainConstraint) or isinstance(other, SetConstraint):
            return False
        elif isinstance(other, InequalityConstraint) or isinstance(other, EqualityConstraint):
            return True
        elif isinstance(other, LessThanConstraint):
            return self.terms < other.terms
        else:
            raise NotImplementedError(f'Cannot compare LessThanConstraint with {type(other)}')


class InclusionConstraint(SetConstraint):
    """
    A FOL-DC inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalVariable', domain_term: 'DomainTerm') -> None:
        """NOTE TODO: For now, this only works with 'SetOfConstants' rather than general 'DomainTerm'
        for the domain term"""
        self.terms: Tuple['LogicalVariable', 'DomainTerm'] = (logical_term, domain_term)

    @property
    def logical_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def domain_term(self) -> 'DomainTerm':
        return self.terms[1]

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        TODO: Include domain substitutions"""
        new_logical_term = substitution[self.logical_term]
        if isinstance(new_logical_term, Constant):
            if isinstance(self.domain_term, SetOfConstants):
                if new_logical_term in self.domain_term:
                    return EmptyConstraint(f'{new_logical_term} in {self.domain_term}')
                else:
                    return FalseConstraint(f'{new_logical_term} not in {self.domain_term}')
            raise ValueError('Cannot yet handle DomainTerm in substitute')
        else:
            logical_var = cast('LogicalVariable', new_logical_term)
            return InclusionConstraint(logical_var, self.domain_term)

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For InclusionConstraint, this means checking if the logical term is a constant and
        the domain term a set of constants without the constant"""
        if isinstance(self.logical_term, Constant) and isinstance(self.domain_term, SetOfConstants):
            return not self.logical_term in self.domain_term
        else:
            return False

    def __eq__(self, other: Any) -> bool:
        """Two inclusion constraints are equal if they mention the same logical and domain terms """
        return isinstance(other, InclusionConstraint) \
               and self.logical_term == other.logical_term \
               and self.domain_term == other.domain_term

    def __hash__(self) -> int:
        """Just using the parent hash function.
        We have to redefine it because we overrode __eq__.

        NOTE: this may cause collisions between InclusionConstraints and NotInclusionConstraints"""
        return super().__hash__()

    def __invert__(self) -> 'SetConstraint':
        """This method overrides the '~' operator.
        I use it to negate the constraint -- e.g. turn = into !=
        NOTE: If this involves a DomainVariable, we return another InclusionConstraint with its complement"""
        if isinstance(self.domain_term, DomainVariable):
            return InclusionConstraint(self.logical_term, self.domain_term.complement)
        else:
            return NotInclusionConstraint(self.logical_term, self.domain_term)

    def __str__(self) -> str:
        element_of_string = ' \u2208 '
        return f'{self.logical_term}{element_of_string}{self.domain_term}'

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        if isinstance(other, DomainConstraint):
            return False
        elif isinstance(other, LogicalConstraint):
            return True
        elif isinstance(other, NotInclusionConstraint):
            return True
        elif isinstance(other, InclusionConstraint):
            return self.terms < other.terms
        else:
            raise NotImplementedError(f'Cannot compare InclusionConstraint with {type(other)}')

class NotInclusionConstraint(SetConstraint):
    """
    A FOL-DC negated inclusion constraint (between a logical term and a domain term).
    This consists of a logical term and a domain term.
    """

    def __init__(self, logical_term: 'LogicalVariable', domain_term: 'DomainTerm') -> None:
        """NOTE: For now, this only works with 'SetOfConstants' rather than general 'DomainTerm'
        for the domain term"""
        self.terms: Tuple['LogicalVariable', 'DomainTerm'] = (logical_term, domain_term)

    @property
    def logical_term(self) -> 'LogicalVariable':
        return self.terms[0]

    @property
    def domain_term(self) -> 'DomainTerm':
        return self.terms[1]

    def substitute(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        TODO: Include domain substitutions"""
        new_logical_term = substitution[self.logical_term]
        if isinstance(new_logical_term, Constant):
            if isinstance(self.domain_term, SetOfConstants):
                if new_logical_term in self.domain_term:
                    return FalseConstraint(f'{new_logical_term} in {self.domain_term}')
                else:
                    return EmptyConstraint(f'{new_logical_term} not in {self.domain_term}')
            raise ValueError('Cannot yet handle DomainTerm in substitute')
        else:
            logical_var = cast('LogicalVariable', new_logical_term)
            return NotInclusionConstraint(logical_var, self.domain_term)

    def contains_contradiction(self) -> bool:
        """Does this constraint contain an obvious contradiction?
        For NotInclusionConstraint, this means checking if the logical term is a constant and the domain term a set of constants containing the constant"""
        if isinstance(self.logical_term, Constant) and isinstance(self.domain_term, SetOfConstants):
            return self.logical_term in self.domain_term
        else:
            return False

    def __eq__(self, other: Any) -> bool:
        """Two not inclusion constraints are equal if they mention the same logical and domain terms """
        return isinstance(other, NotInclusionConstraint) \
               and self.logical_term == other.logical_term \
               and self.domain_term == other.domain_term

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

    def __lt__(self, other: Any) -> bool:
        if isinstance(other, DomainConstraint):
            return False
        elif isinstance(other, LogicalConstraint):
            return True
        elif isinstance(other, InclusionConstraint):
            return False
        elif isinstance(other, NotInclusionConstraint):
            return self.terms < other.terms
        else:
            raise NotImplementedError(f'Cannot compare NotInclusionConstraint with {type(other)}')


class Substitution:
    """A FOL substitution.
    This contains a dictionary of logical variables and their substitutions (which are terms)
    TODO: expand this to work with domain variables
    """
    def __init__(self, variable_term_pairs: Iterable[VarTermPair]) -> None:
        """The substitution dict is private and shouldn't be changed after creation"""
        self._substitution_dict = {var: term for var, term in sorted(variable_term_pairs)}

    def __getitem__(self, query: 'LogicalTerm') -> 'LogicalTerm':
        """Return the term associated with a given (logical) variable
        NOTE: returns the query if it is not present in the substitution"""
        if not isinstance(query, LogicalVariable):
            return query
        value = self._substitution_dict.get(query)
        if value is None:
            value = query
        return value

    def mappings(self) -> Iterable[VarTermPair]:
        """Return an iterator of the mappings in this substitution."""
        return self._substitution_dict.items()

    def to_constraint_set(self) -> 'ConstraintSet':
        """Convert a substitution to an equivalent constraint set.
        This involves simply defining a constraint for each mapping."""
        constraints: List['Constraint'] = []
        for mapping in self:
            if isinstance(mapping[1], LogicalVariable):
                constraint = EqualityConstraint(mapping[0], mapping[1])
            elif isinstance(mapping[1], Constant):
                constraint = InclusionConstraint(mapping[0], SetOfConstants([mapping[1]]))
            constraints.append(constraint)
        return ConstraintSet(constraints)

    def __iter__(self):
        return iter(self.mappings())

    def __eq__(self, other: Any) -> bool:
        """Two substitutions are equal if they have all the same (variable, term) pairs"""
        return isinstance(other, Substitution) and self._substitution_dict == other._substitution_dict

    def __hash__(self) -> int:
        """Hashing so I can remove duplicates with sets"""
        return hash(tuple((key, val) for key, val in self._substitution_dict.items()))

    def __str__(self) -> str:
        rightarrow_string = ' \u2192 '
        substitution_strings = [str(v) + rightarrow_string + str(t) for v, t in self._substitution_dict.items()]
        return f"{{{', '.join(substitution_strings)}}}"

    def __repr__(self) -> str:
        return self.__str__()



class EquivalenceClass:
    """A class implementing an equivalence class between logical terms in FOL-DC.
    """
    def __init__(self, members: Iterable['LogicalTerm']) -> None:
        """Instantiate the private set of members that belong to this equivalence class"""
        self._members = frozenset(members)
        self._variables = frozenset(t for t in members if isinstance(t, LogicalVariable))
        self._constants = frozenset(t for t in members if isinstance(t, Constant))

    @property
    def members(self) -> FrozenSet['LogicalTerm']:
        return self._members

    @property
    def variables(self) -> FrozenSet['LogicalVariable']:
        return self._variables

    @property
    def constants(self) -> FrozenSet['Constant']:
        return self._constants

    def overlaps(self, other: 'EquivalenceClass') -> bool:
        """Returns True if there is an element shared by this EquivalenceClass and the other"""
        return len(self.members.intersection(other.members)) > 0

    def is_consistent(self) -> bool:
        """If this equivalence class contains an obvious contradiction (i.e. two different constants)
        then it is inconsistent and we return False"""
        for term1 in self.members:
            if isinstance(term1, Constant):
                for term2 in self.members:
                    if isinstance(term2, Constant) and term1 != term2:
                        return False
        return True

    def is_inconsistent(self) -> bool:
        """Negation of is_consistent (for convenience)"""
        return not self.is_consistent

    def join(self, other: 'EquivalenceClass') -> 'EquivalenceClass':
        """Create a larger equivalence class by adding the members of another"""
        new_members = self.members.union(other.members)
        return EquivalenceClass(new_members)

    def get_variables_only(self) -> 'VariableEquivalenceClass':
        """Take an equivalence class and remove all constants, returning a subclass"""
        return VariableEquivalenceClass(self._variables)


    def get_shared_domain_from_cs(self, cs: 'ConstraintSet') -> 'ProperDomain':
        """Get the shared domain for this equivalence class of variables according to a particular constraint set
        NOTE: This now works only with ProperDomains
        TODO TODO: Add excluded domain part back in"""
        included_domains = []
        excluded_domains = []
        for constraint in cs:
            if isinstance(constraint, InclusionConstraint) and constraint.logical_term in self:
                if isinstance(constraint.domain_term, ProperDomain):
                    included_domains.append(constraint.domain_term)
            elif isinstance(constraint, NotInclusionConstraint) and constraint.logical_term in self:
                if isinstance(constraint.domain_term, ProperDomain):
                    excluded_domains.append(constraint.domain_term)
        shared_domain = ProperDomain.intersect_all(*included_domains)
        # check that none of the excluded domains are supersets of the shared domain
        for excluded_domain in excluded_domains:
            if excluded_domain.is_strict_superset_of(shared_domain) or excluded_domain == shared_domain:
                return EmptyDomain(f'Excluded {excluded_domain} is superset of shared {shared_domain}')
            elif excluded_domain.is_strict_subset_of(shared_domain):
                # NOTE TODO: trying out a complement domain thing
                if excluded_domain.parent_domain == shared_domain and isinstance(excluded_domain, DomainVariable):
                    return excluded_domain.complement
                raise ValueError(f'Complex shared domain! Shared {shared_domain}, excluded {excluded_domain}')
        return shared_domain

    def is_root_in_cnf(self, cnf: 'CNF') -> bool:
        """Determine whether this equivalence class is root for a CNF
        (i.e. each variable appears in every literal of its clause)"""
        # the variables that appear in a clause must all be root
        for clause in cnf.clauses:
            if self.members.intersection(clause.literal_variables) != self.members.intersection(clause.root_variables):
                return False
        return True

    def _get_constant_or_random_variable(self) -> 'LogicalTerm':
        """If there's a constant in this equivalence class, return that.
        Otherwise, return an arbitrary variable"""
        for element in sorted(self.members):
            if isinstance(element, Constant):
                return element
        return element

    def _partition_overlapping_disjoint_classes(self,
                                                other_eq_classes: Set['TEC']
                                                ) -> Tuple[Set['TEC'], Set['TEC']]:
        """Return two lists of equivalence classes, one for those classes that overlap with this
        class, and another for the rest"""
        overlapping: Set['TEC']
        disjoint: Set['TEC']
        overlapping, disjoint = set(), set()
        for eq_class in other_eq_classes:
            if eq_class.overlaps(self):
                overlapping.add(eq_class)
            else:
                disjoint.add(eq_class)
        return overlapping, disjoint


    def __iter__(self):
        return iter(self.members)

    def __str__(self) -> str:
        return str(self._members)

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """Equivalence classes can be compared by sorting their members"""
        if isinstance(other, EquivalenceClass):
            return sorted(self.members) < sorted(other.members)
        else:
            raise NotImplementedError(f'Cannot compare EquivalenceClass and {type(other)}')


class VariableEquivalenceClass(EquivalenceClass):
    """A subclass of EquivalenceClass that can contain only logical variables"""
    def __init__(self, members: Iterable['LogicalVariable']) -> None:
        self._members = frozenset(members)

class EquivalenceClasses(Generic[TEC]):
    """This is wrapper class around a set of 'EquivalenceClass' objects.
    NOTE: It uses generic typing to work with either EquivalenceClass or VariableEquivalenceClass"""
    def __init__(self, eq_classes: Iterable['TEC']) -> None:
        self._eq_classes = set(eq_classes)

    @property
    def classes(self) -> Set[TEC]:
        return self._eq_classes

    def make_self_consistent(self) -> 'EquivalenceClasses[TEC]':
        """Take an initial set of equivalence classes and iterate them until they are
        self-consistent, i.e. no overlapping equivalence classes.
        NOTE: Could be a little more efficient if we checked for inconsistency of constants 
        earlier, but that would make it less clean.
        TODO: Decide if this makes sense as the place for this function"""
        remaining_eq_classes = self.classes.copy()
        final_eq_classes: Set[TEC] = set()
        while len(remaining_eq_classes) > 0:
            current_eq_class = remaining_eq_classes.pop()

            current_class_changed = True
            while current_class_changed:
                # split the equivalence classes into those that overlap with the current class and those that don't
                overlapping, disjoint = current_eq_class._partition_overlapping_disjoint_classes(remaining_eq_classes)
                # merge the overlapping eq classes into one big eq class
                merger = lambda eq_class1, eq_class2: eq_class1.join(eq_class2)
                current_eq_class = reduce(merger, overlapping, current_eq_class)
                remaining_eq_classes = disjoint
                current_class_changed = (len(overlapping) > 0)

            final_eq_classes.add(current_eq_class)
        # return as an instance of TEC (either 
        return self.__class__(final_eq_classes)

    def to_substitution(self) -> 'Substitution':
        """Convert a list of equivalence classes into a Substitution of variables that conveys
        the same information."""
        var_term_pairs: List[VarTermPair] = []
        # we process each equivalence class in turn, generating var-term pairs
        for eq_class in sorted(self.classes):
            # choose what to map everything to
            mapping_target = eq_class._get_constant_or_random_variable()
            # now we get term-var pairs for all terms that are NOT the target
            for term in sorted(eq_class.members):
                if term != mapping_target:
                    if isinstance(term, Constant):
                        raise ValueError('There should have been at most 1 constant')
                    var = cast('LogicalVariable', term) # hack for type-checking
                    var_term_pairs.append((var, mapping_target))
        return Substitution(var_term_pairs)

    def to_constraint_set(self) -> 'ConstraintSet':
        """Convert a sequence of equivalence classes into a substitution.
        We go via a subtitution because we already have functions for that"""
        substitution = self.to_substitution()
        cs = substitution.to_constraint_set()
        return cs

    def consistent_with_variable_inequality_constraints(self, cs: 'ConstraintSet') -> bool:
        """If the equality classes conflict with the inequality constraints between variables, returns False"""
        for eq_class in self.classes:
            for constraint in cs.inequality_constraints:
                if constraint.left_term in eq_class and constraint.right_term in eq_class: 
                    return False
        return True

    def consistent_with_constant_inequality_constraints(self, cs: 'ConstraintSet') -> bool:
        """If the equality classes conflict with inequalities between constants, return False"""
        for eq_class in self.classes:
            equal_constants: Set['Constant'] = set()
            unequal_constants: Set['Constant'] = set()
            for ic in cs.inclusion_constraints:
                if isinstance(ic.domain_term, SetOfConstants) and ic.domain_term.size() == 1 and ic.logical_term in eq_class:
                    equal_constants.add(get_element_of_set(ic.domain_term.constants))
            for nc in cs.notinclusion_constraints:
                if isinstance(nc.domain_term, SetOfConstants) and nc.domain_term.size() == 1 and nc.logical_term in eq_class:
                    unequal_constants.add(get_element_of_set(nc.domain_term.constants))

            if len(equal_constants) > 1:
                return False  # can't be equal to more than one constant
            if len(equal_constants.intersection(unequal_constants)) >= 1:
                return False  # can't be equal to constants we have an inequality with
        return True

    def consistent_with_set_constraints(self, cs: 'ConstraintSet') -> bool:
        """If any equality class WITH A PROPER DOMAIN has an empty domain, then it is not consistent (inclusion constraint condition is so free variables don't mess it up)
        SetOfConstants are now handled in separately"""
        if len(cs.inclusion_constraints) == 0:
            return True  # no inclusion constraints to conflict with
        for eq_class in self.classes:
            variables_with_domains = set(c.logical_term for c in cs.inclusion_constraints if c.logical_term in eq_class and isinstance(c.domain_term, ProperDomain))

            if len(variables_with_domains) == 0:
                return True # trivially consistent

            reduced_eq_class = EquivalenceClass(variables_with_domains)
            shared_domain = reduced_eq_class.get_shared_domain_from_cs(cs)
            if shared_domain.size() == 0:
                return False
        return True

    def __iter__(self):
        return iter(self.classes)

    def __str__(self) -> str:
        return f'ECs{str(self.classes)}'

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """EquivalenceClasses can be compared by sorting their member classes"""
        if isinstance(other, EquivalenceClasses):
            return sorted(self.classes) < sorted(other.classes)
        else:
            raise NotImplementedError(f'Cannot compare EquivalenceClasses and {type(other)}')
