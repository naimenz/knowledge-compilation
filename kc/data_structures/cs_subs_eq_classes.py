"""
Classes for domain constraints and constraint sets AND NOW Substitutions, because they depend on each other. 
AND NOW ALSO EQUIVALENCE CLASSES
NOTE: at the moment only the types of constraints used in Chapter 4 of the PhD
are considered.
"""

from kc.data_structures import LogicalVariable, SetOfConstants, Constant, DomainTerm

from abc import ABC, abstractmethod
from functools import reduce

from typing import List, Iterable, Any, FrozenSet, Tuple, Union, Optional, Sequence, Set
from typing import TypeVar, cast, Generic
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import CNF, LogicalVariable, LogicalTerm

VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']
TEC = TypeVar('TEC', bound='EquivalenceClass') # this is a type var so I can work with eq classes or var eq classes

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

    def is_non_empty(self) -> bool:
        return len(self.constraints) > 0

    def apply_substitution(self, substitution: 'Substitution') -> 'ConstraintSet':
        """Create a ConstraintSet by applying a Substitution to the constraints
        NOTE: this may be counter-intuitive because applying a substitution
        to an EqualityConstraint may produce an InclusionConstraint"""
        new_constraints: List['Constraint'] = []
        for constraint in self.constraints:
            new_constraint = constraint.apply_substitution(substitution)
            new_constraints.append(new_constraint)

        if FalseConstraint() in new_constraints:
            raise ValueError(f'Substitution {substitution} made cs unsatisfiable!')

        # remove trivial constraints (always true)
        filtered_constraints = [c for c in new_constraints if c != EmptyConstraint()]

        return ConstraintSet(filtered_constraints)

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
        if not eq_classes.consistent_with_inequality_constraints(self):
            return False
        # Then we construct the possible domains for each equivalence class and
        # check that they are non-empty
        if not eq_classes.consistent_with_set_constraints(self):
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

    @abstractmethod 
    def apply_substitution(self, substitution: 'Substitution') -> 'Constraint':
        """We need to be able to apply substitutions to constraints"""
        pass

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

class EmptyConstraint(Constraint):
    """This is a special class for a constraint that is trivially satisfied (i.e. always true)"""
    def apply_substitution(self, substitution: 'Substitution') -> 'EmptyConstraint':
        """Substitution does nothing because no terms"""
        return EmptyConstraint()
    
    def contains_contradiction(self) -> bool:
        """Never contains a contradiction because it's always true"""
        return False

    def __hash__(self) -> int:
        """All EmptyConstraints are the same"""
        return hash('EmptyConstraint')


class FalseConstraint(Constraint):
    """This is a special class for a constraint that is trivially UNsatisfied (i.e. always false)"""
    def apply_substitution(self, substitution: 'Substitution') -> 'FalseConstraint':
        """Substitution does nothing because no terms"""
        return FalseConstraint()
    
    def contains_contradiction(self) -> bool:
        """Always contains a contradiction because it's always false"""
        return True

    def __hash__(self) -> int:
        """All EmptyConstraints are the same"""
        return hash('FalseConstraint')


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

    def apply_substitution(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        NOTE: This will return an InclusionConstraint if ONE of the new terms is a constant.
        TODO: Figure out how to abstract this up the the base class"""
        new_left_term = substitution[self.left_term]
        new_right_term = substitution[self.right_term]
        if isinstance(new_left_term, Constant):
            if isinstance(new_right_term, Constant):
                if new_left_term == new_right_term:
                    return EmptyConstraint()
                else:
                    return FalseConstraint()
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

    def apply_substitution(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        NOTE: This will return a NotInclusionConstraint if ONE of the new terms is a constant.
        """
        new_left_term = substitution[self.left_term]
        new_right_term = substitution[self.right_term]
        if isinstance(new_left_term, Constant):
            if isinstance(new_right_term, Constant):
                if new_left_term != new_right_term:
                    return EmptyConstraint()
                else:
                    return FalseConstraint()
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

    def apply_substitution(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        TODO: Include domain substitutions"""
        new_logical_term = substitution[self.logical_term]
        if isinstance(new_logical_term, Constant):
            if isinstance(self.domain_term, SetOfConstants):
                if new_logical_term in self.domain_term:
                    return EmptyConstraint()
                else:
                    return FalseConstraint()
        else:
            logical_var = cast('LogicalVariable', new_logical_term)
            return InclusionConstraint(logical_var, self.domain_term)

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

    def apply_substitution(self, substitution: 'Substitution') -> 'Constraint':
        """Apply a substitution to the constraint, returning a new constraint.
        TODO: Include domain substitutions"""
        new_logical_term = substitution[self.logical_term]
        if isinstance(new_logical_term, Constant):
            if isinstance(self.domain_term, SetOfConstants):
                if new_logical_term in self.domain_term:
                    return FalseConstraint()
                else:
                    return EmptyConstraint()
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



class Substitution:
    """A FOL substitution.
    This contains a dictionary of logical variables and their substitutions (which are terms)
    TODO: expand this to work with domain variables
    """
    def __init__(self, variable_term_pairs: Iterable[VarTermPair]) -> None:
        """The substitution dict is private and shouldn't be changed after creation"""
        self._substitution_dict = {var: term for var, term in variable_term_pairs}

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
            constraint = EqualityConstraint(mapping[0], mapping[1])
            constraints.append(constraint)
        return ConstraintSet(constraints)

    def __iter__(self):
        return iter(self.mappings())

    def __eq__(self, other: Any) -> bool:
        """Two substitutions are equal if they have all the same (variable, term) pairs"""
        if not isinstance(other, Substitution):
            return False
        return self._substitution_dict == other._substitution_dict

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
        self._constants = frozenset(t for t in members if isinstance(t, LogicalVariable))

    @property
    def members(self) -> FrozenSet['LogicalTerm']:
        return self._members

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
        """Take an equivalence class and remove all constants"""
        variables = [term for term in self.members if isinstance(term, LogicalVariable)]
        return VariableEquivalenceClass(variables)


    def get_shared_domain_from_cs(self, cs: 'ConstraintSet') -> 'SetOfConstants':
        """Get the shared domain for this equivalence class of variables according to a particular constraint set
        NOTE: Because we don't have inequality constraints involving constants, this covers everything"""
        included_domains = []
        excluded_domains = []
        for constraint in cs:
            if isinstance(constraint, InclusionConstraint) and constraint.logical_term in self:
                included_domains.append(constraint.domain_term)
            elif isinstance(constraint, NotInclusionConstraint) and constraint.logical_term in self:
                excluded_domains.append(constraint.domain_term)
        # hack for type checking for now since we only work with SetOfConstants
        included_domain = cast('SetOfConstants', DomainTerm.intersection(*included_domains)) 
        excluded_domain = cast('SetOfConstants', DomainTerm.union(*excluded_domains)) 
        shared_domain = included_domain.difference(excluded_domain)
        return shared_domain

    def is_root_in_cnf(self, cnf: 'CNF') -> bool:
        """Determine whether this equivalence class is root for a CNF
        (i.e. each variable appears in every literal of its clause)"""
        # the variables that appear in a clause must all be root
        for clause in cnf.clauses:
            if self.members.intersection(clause.all_literal_variables) != self.members.intersection(clause.root_variables):
                return False
        return True

    def _get_constant_or_random_variable(self) -> 'LogicalTerm':
        """If there's a constant in this equivalence class, return that.
        Otherwise, return an arbitrary variable"""
        for element in self.members:
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


class VariableEquivalenceClass(EquivalenceClass):
    """A subclass of EquivalenceClass that can contain only logical variables"""
    def __init__(self, members: Iterable['LogicalVariable']) -> None:
        super().__init__(members)


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
        for eq_class in self.classes:
            # choose what to map everything to
            mapping_target = eq_class._get_constant_or_random_variable()
            # now we get term-var pairs for all terms that are NOT the target
            for term in eq_class.members:
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

    def consistent_with_inequality_constraints(self, cs: 'ConstraintSet') -> bool:
        """If the equality classes conflict with the inequality constraints, returns False
        NOTE: assuming the inequality constraints only contain logical variables"""
        for eq_class in self.classes:
            for constraint in cs:
                if isinstance(constraint, InequalityConstraint):
                    if constraint.left_term in eq_class and constraint.right_term in eq_class: 
                        return False
        return True

    def consistent_with_set_constraints(self, cs: 'ConstraintSet') -> bool:
        """If any equality class has an empty domain, then it is not consistent
        NOTE: For now, only works with SetOfConstants, not DomainVariable"""
        for eq_class in self.classes:
            shared_domain = eq_class.get_shared_domain_from_cs(cs)
            if shared_domain.size == 0:
                return False
        return True

    def __iter__(self):
        return iter(self.classes)

    def __str__(self) -> str:
        return f'ECs{str(self.classes)}'

    def __repr__(self) -> str:
        return self.__str__()
