"""
Classes for clauses in FOL-DC.
This includes constrained AND unconstrained clauses.

TODO: Figure out if the inheritance structure for UnitClause and ConstrainedAtom makes sense.
"""

from kc.data_structures import Literal, Atom, LogicalVariable, Constant, ConstraintSet, InequalityConstraint, NotInclusionConstraint, SetOfConstants, CNF

from functools import reduce
from abc import ABC, abstractmethod

from typing import List, TypeVar, Iterable, Any, Sequence, Set, Optional
from typing import cast 
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import Atom, Substitution, ConstraintSet, EquivalenceClasses, Constraint

# Type variable for arbitrary clauses so I can reuse substitute
C = TypeVar('C', bound='ConstrainedClause') 

class Clause(ABC):
    """Abstract base class for constrained and unconstrained clauses"""
    def __init__(self, literals: Iterable['Literal']) -> None:
        self.literals = frozenset(literals)

    @abstractmethod
    def propagate_equality_constraints(self: 'Clause') -> 'Clause':
        """Propagate the equality constraints through the clause by building a substitution
        from them and applying it until convergence.
        NOTE: this is only important for ConstrainedClauses, but needs to be callable on both"""
        pass

    @abstractmethod
    def substitute(self: 'Clause', substitution: 'Substitution') -> 'Clause':
        pass

    @abstractmethod
    def get_constrained_atoms(self) -> List['ConstrainedAtom']:
        """Even though UnconstrainedClauses don't have constraints,
        they still need to produce constrained atoms (with empty bound_vars and
        constraints)"""
        pass

    @abstractmethod
    def get_constrained_literals(self) -> List['UnitClause']:
        """Same as get_constrained_atoms but for literals.  """
        pass

    @abstractmethod
    def to_unit_clause(self) -> 'UnitClause':
        """Convert this clause to a UnitClause object"""
        pass

    def is_independent_from_other_clause(self, other_clause: 'Clause') -> bool:
        """Is this clause independent of the other clause?"""
        if self.has_no_literals() or other_clause.has_no_literals():
            return True

        for c_atom in self.get_constrained_atoms():
            for other_c_atom in other_clause.get_constrained_atoms():
                if c_atom.constrained_atoms_unify(other_c_atom):
                    return False
        return True

    def is_subsumed_by_literal(self, subsumer: 'UnitClause') -> bool:
        """Returns true if this Clause is subsumed by the constrained literal.
        NOTE: We now say that EVERY literal in self has to be subsumed by some literal in the subsumer
        TODO: check this is right"""
        if self.has_no_literals():
            return True

        # every literal in the subsumer must subsume some literal in self
        for c_literal in self.get_constrained_literals():
            if c_literal.is_subsumed_by_literal(subsumer):
                return True
        return False

    def has_no_literals(self) -> bool:
        """This is called "isConditionalContradiction" in Forclift.
        I think this means that the clause contains no literals and so
        its grounding is empty, meaning it is independent of everything."""
        return len(self.literals) == 0

    @abstractmethod
    def is_tautology(self) -> bool:
        pass

    def is_contradiction(self) -> bool:
        """Is this clause always false? For now, just check if
        its constraint set is satisfiable and it contains a literal
        and its negation
        NOTE: I don't think that a clause can be a contradiction because
        it doesn't contain ANDs
        TODO: check whether clauses can be contradictions"""
        return False

    @property
    def literal_variables(self) -> Set['LogicalVariable']:
        """Return all variables that appear in any literal of the clause"""
        all_vars: Set['LogicalVariable'] = set()
        for literal in self.literals:
            all_vars = all_vars.union(literal.variables)
        return all_vars

    @property
    @abstractmethod
    def all_variables(self) -> Set['LogicalVariable']:
        """Return ALL variables that appear in the clause (this includes bound vars
        and cs for ConstrainedClauses)"""

    @property
    def root_variables(self) -> Set['LogicalVariable']:
        """Return the root variables of this clause 
        (i.e. the variables that appear in every literal of the clause)"""
        # first, collect ALL variables that appear in any literal
        root_vars = self.literal_variables
        # then get only the ones that appear in every literal
        for literal in self.literals:
            root_vars = root_vars.intersection(literal.variables)
        return root_vars

    @property
    def constants(self) -> Set['Constant']:
        """Extract just the constants as logical terms from the clause.
        NOTE: Since we do not allow constants in LogicalConstraints, 
        we only have to look at the literals"""
        all_constants: Set['Constant'] = set()
        for literal in self.literals:
            all_constants = all_constants.union(literal.constants)
        return all_constants


class UnconstrainedClause(Clause):
    """An FOL unconstrained clause.
    This consists of a set of FOL literals, which form a disjunction
    """

    def get_constrained_atoms(self) -> List['ConstrainedAtom']:
        """Even though UnconstrainedClauses don't have constraints,
        they still need to produce constrained atoms (with empty bound_vars and
        constraints)"""
        constrained_atoms = []
        for literal in self.literals:
            empty_bound_vars: Set['LogicalVariable'] = set()
            empty_cs = ConstraintSet([])
            positive_literal = Literal(literal.atom, True)
            constrained_atom = ConstrainedAtom([positive_literal], empty_bound_vars, empty_cs)
            constrained_atoms.append(constrained_atom)
        return constrained_atoms

    def get_constrained_literals(self) -> List['UnitClause']:
        """Even though UnconstrainedClauses don't have constraints,
        they still need to produce constrained literals (with empty bound_vars and
        constraints)"""
        constrained_literals = []
        for literal in self.literals:
            empty_bound_vars: Set['LogicalVariable'] = set()
            empty_cs = ConstraintSet([])
            constrained_literal = UnitClause([literal], empty_bound_vars, empty_cs)
            constrained_literals.append(constrained_literal)
        return constrained_literals

    def is_tautology(self) -> bool:
        """Is this clause always true? 
        For now, just check if it contains a literal and its negation"""
        for literal in self.literals:
            if ~literal in self.literals:
                return True
        return False

    @property
    def all_variables(self) -> Set['LogicalVariable']:
        """Return ALL variables. For an UnconstrainedClause, this is just
        the literal variables"""
        return self.literal_variables

    def to_unit_clause(self) -> 'UnitClause':
        """Convert this clause to a UnitClause object.
        NOTE: throws an error if this is not possible"""
        if len(self.literals) != 1:
            raise ValueError('Cannot convert UnconstrainedClause with {len(self.literals)} literals to a UnitClause')
        empty_bound_vars: Set['LogicalVariable'] = set()
        empty_cs = ConstraintSet([])
        return UnitClause(self.literals, empty_bound_vars, empty_cs)

    def propagate_equality_constraints(self: 'UnconstrainedClause') -> 'UnconstrainedClause':
        """Propagate the equality constraints through the clause by building a substitution
        from them and applying it until convergence."""
        return self

    def substitute(self: 'UnconstrainedClause', substitution: 'Substitution') -> 'UnconstrainedClause':
        """Return a new UnconstrainedClause, the result of applying substitution to this UnconstrainedClause"""
        return self.__class__(literal.substitute(substitution) for literal in self.literals)

    def __eq__(self, other: Any) -> bool:
        """Two unconstrained clauses are equal if they have the same literals"""
        return isinstance(other, UnconstrainedClause) and self.literals == other.literals

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
            literals: Iterable['Literal'],
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        self.literals = frozenset(literals)
        self.bound_vars = frozenset(bound_vars)
        self.cs = cs

    def substitute(self: 'C', substitution: 'Substitution') -> 'C':
        """Return a new ConstrainedClause, the result of applying substitution to this ConstrainedClause
        NOTE: For now we allow substitution of constants to bound vars by just having one fewer bound var"""
        new_literals = [literal.substitute(substitution) for literal in self.literals]
        new_cs = self.cs.substitute(substitution)
        _new_bound_vars = [substitution[var] for var in self.bound_vars if isinstance(substitution[var], LogicalVariable)]
        ## allowing substitution of constants to bound vars at the moment
        # assert(all(isinstance(term, LogicalVariable) for term in _new_bound_vars)) 
        new_bound_vars = cast(List['LogicalVariable'], _new_bound_vars) # hack for type checking
        return self.__class__(new_literals, new_bound_vars, new_cs)

    def propagate_equality_constraints(self: 'C') -> 'C':
        """Propagate the equality constraints through the clause by building a substitution
        from them and applying it until convergence."""
        var_eq_classes = self.cs.get_var_eq_classes()
        sub = var_eq_classes.to_substitution()
        # repeatedly apply the substitution until convergence
        old_clause = self
        while True:
            new_clause = old_clause.substitute(sub)
            if new_clause != old_clause:
                old_clause = new_clause
            else:
                return new_clause

    @property
    def all_variables(self) -> Set['LogicalVariable']:
        """Return ALL variables, regardless of whether they are bound, in the cs,
        or in the literals"""
        literal_vars = self.literal_variables
        bound_vars = self.bound_vars
        constraint_vars = self.constraint_variables 
        return literal_vars.union(bound_vars, constraint_vars)


    @property
    def constraint_variables(self) -> Set['LogicalVariable']:
        """Extract just the variables from each constraint in the constraint set
        NOTE: I had to duplicate this here to avoid circular imports"""
        logical_variables: Set['LogicalVariable'] = set()
        for constraint in self.cs:
            for term in constraint.terms:
                if isinstance(term, LogicalVariable):
                    logical_variables.add(term)
        return logical_variables

    def get_free_variables(self) -> Set['LogicalVariable']:
        """Extract just the free variables from this clause"""
        free_variables: Set['LogicalVariable'] = set()
        for variable in self.constraint_variables.union(self.literal_variables):
            if variable not in self.bound_vars:
                free_variables.add(variable)
        return free_variables

    # TODO: THINK ABOUT FREE VARIABLE CASE
    def get_constant_or_free_inequalities(self) -> Set['Constraint']:
        """Get inequalities that are between a bound variable and a constant OR FREE VARIABLE in this clause.
        NOTE: these will be NotInclusionConstraints because of how I've implemented those
        NOTE 2: For now, just treating free variables exactly like constnats"""
        ineq_constraints: Set['Constraint'] = set()
        # get constant inequalities
        for set_constraint in self.cs.set_constraints:
            if isinstance(set_constraint, NotInclusionConstraint):
                domain_term = set_constraint.domain_term
                if isinstance(domain_term, SetOfConstants) and domain_term.size() == 1 and set_constraint.logical_term in self.bound_vars:
                    ineq_constraints.add(set_constraint)
        # get free_variable inequalities
        free_variables = self.get_free_variables()
        bound_variables = self.bound_vars
        for logical_constraint in self.cs.logical_constraints:
            if isinstance(logical_constraint, InequalityConstraint):
                if (logical_constraint.left_term in free_variables and logical_constraint.right_term in bound_variables) \
                or (logical_constraint.right_term in free_variables and logical_constraint.left_term in bound_variables):
                    ineq_constraints.add(logical_constraint)
        return ineq_constraints

    def get_bound_variable_inequalities(self) -> Set['InequalityConstraint']:
        """Get inequalities that are between bound variables in this clause"""
        ineq_constraints = set()
        for constraint in self.cs.logical_constraints:
            if isinstance(constraint, InequalityConstraint):
                if constraint.left_term in self.bound_vars and constraint.right_term in self.bound_vars:
                    ineq_constraints.add(constraint)
        return ineq_constraints

    def is_tautology(self) -> bool:
        """Is this clause always true? For now, just check if
        its constrainst set is satisfiable and it contains a literal
        and its negation"""
        if not self.cs.is_satisfiable():
            return True

        for literal in self.literals:
            if ~literal in self.literals:
                return True
        return False

    def is_contradiction(self) -> bool:
        """Is this clause always false? For now, just check if
        its constraint set is satisfiable and it contains a literal
        and its negation
        NOTE: I don't think that a clause can be a contradiction because
        it doesn't contain ANDs
        TODO: check whether clauses can be contradictions"""
        return False

    def has_no_literals(self) -> bool:
        """This is called "isConditionalContradiction" in Forclift.
        I think this means that the clause contains no literals and so
        its grounding is empty, meaning it is independent of everything."""
        return len(self.literals) == 0

    def get_constrained_atoms(self) -> List['ConstrainedAtom']:
        """For this clause, return a list of all the constrained atoms in the clause"""
        constrained_atoms = []
        for literal in self.literals:
            constrained_atom = self._build_constrained_atom_from_literal(literal)
            constrained_atoms.append(constrained_atom)
        return constrained_atoms

    def get_constrained_literals(self) -> List['UnitClause']:
        """For this clause, return a list of all the constrained atoms in the clause"""
        constrained_literals = []
        for literal in self.literals:
            constrained_literal = self._build_constrained_literal_from_literal(literal)
            constrained_literals.append(constrained_literal)
        return constrained_literals

    def _build_constrained_atom_from_literal(self, literal: 'Literal') -> 'ConstrainedAtom':
        """Build a constrained atom given a specific literal using the bound variables
        and constraints of this clause"""
        atom = literal.atom 
        positive_literal = Literal(atom, True)
        # constrained atoms subclass ConstrainedClause, so need to be built of an
        # unconstrained clause
        constrained_atom = ConstrainedAtom([positive_literal], self.bound_vars, self.cs)
        return constrained_atom

    def _build_constrained_literal_from_literal(self, literal: 'Literal') -> 'UnitClause':
        """Build a constrained atom given a specific literal using the bound variables
        and constraints of this clause"""
        # constrained atoms subclass ConstrainedClause, so need to be built of an
        # unconstrained clause
        constrained_literal = UnitClause([literal], self.bound_vars, self.cs)
        return constrained_literal

    def to_unit_clause(self) -> 'UnitClause':
        """Convert this clause to a UnitClause object.
        NOTE: throws an error if this is not possible"""
        if len(self.literals) != 1:
            raise ValueError('Cannot convert UnconstrainedClause with {len(self.literals)} literals to a UnitClause')
        return UnitClause(self.literals, self.bound_vars, self.cs)

    def __eq__(self, other: Any) -> bool:
        """Two constrained literals are equal if they have the same unconstrained literals, the same constraint sets,
         and the same bound variables"""
        return isinstance(other, ConstrainedClause) \
               and self.literals == other.literals \
               and self.bound_vars == other.bound_vars \
               and self.cs == other.cs

    def __hash__(self) -> int:
       return hash((self.literals, self.bound_vars, self.cs))

    def __str__(self) -> str:
        bound_vars_strs = [str(var) for var in self.bound_vars]
        for_all_string = '\u2200'
        return f"{for_all_string}{{{', '.join(bound_vars_strs)}}}, {self.cs} : {UnconstrainedClause(self.literals)}"

    def __repr__(self) -> str:
        return self.__str__()


class UnitClause(ConstrainedClause):
    """AN FOL-DC unit clause.
    This is a constrained clause with a single literal
    """
    def __init__(self,
            literals: Iterable['Literal'],
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(tuple(literals)) == 1) # ensure that this is a unit clause
        self.literal = list(literals)[0] # convert to 1-item list and get item
        super(UnitClause, self).__init__(literals, bound_vars, cs)

    def to_c_atom(self) -> 'ConstrainedAtom':
        """Return this as a constrained atom (throw away polarity)"""
        return ConstrainedAtom([Literal(self.literal.atom, True)], self.bound_vars, self.cs)

    def is_subsumed_by_literal(self, subsumer: 'UnitClause') -> bool:
        """Check subsumption of literals.
        This is the same as for atoms but with an additional check of their polarities"""
        if self.literal.polarity == subsumer.literal.polarity \
        and subsumer.to_c_atom().subsumes(self.to_c_atom()):
            return True
        else:
            return False


class ConstrainedAtom(UnitClause):
    """An FOL-DC constrained atom.
    This is a constrained clause with a single atom (positive literal).
    """
    def __init__(self,
            literals: Iterable['Literal'],
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        super(ConstrainedAtom, self).__init__(literals, bound_vars, cs)
        assert(self.literal.polarity) # ensure that it is not negated

    @property
    def atom(self) -> 'Atom':
        """Get the atom from the unconstrained clause"""
        return self.literal.atom


    def constrained_atoms_unify(self: 'ConstrainedAtom', other_c_atom: 'ConstrainedAtom') -> bool:
        """Returns True if there is a substitution that unifies this constrained atom (self) and other_c_atom, otherwise False."""
        return not self.get_constrained_atom_mgu_eq_classes(other_c_atom) is None

    def get_constrained_atom_mgu_eq_classes(self: 'ConstrainedAtom',
                                              other_c_atom: 'ConstrainedAtom'
                                              ) -> Optional['EquivalenceClasses']:
        """Get the mgu of this constrained atom with another.
        This is the same as the mgu for two unconstrained atoms, with an additional check 
        to see if the constraint sets conjoined with the mgu are satisfiable."""
        unconstrained_mgu = self.atom.get_unconstrained_atom_mgu_eq_classes(other_c_atom.atom)
        if unconstrained_mgu is None:
            return None
        cs_mgu = unconstrained_mgu.to_constraint_set()
        combined_constraint_set = self.cs.join(other_c_atom.cs).join(cs_mgu)
        # print(f"DEBUG {combined_constraint_set=}")
        if combined_constraint_set.is_satisfiable():
            return unconstrained_mgu
        else:
            return None

    def get_constrained_atom_mgu_substitution(self: 'ConstrainedAtom',
                                              other_c_atom: 'ConstrainedAtom'
                                              ) -> Optional['Substitution']:
        """Get the mgu of two constrained atoms.
        This is the same as the mgu for two unconstrained atoms, with an additional check 
        to see if the constraint sets conjoined with the mgu are satisfiable."""
        eq_classes = self.get_constrained_atom_mgu_eq_classes(other_c_atom)
        if not eq_classes is None:
            return eq_classes.to_substitution()
        else:
            return None

    def independent_or_subsumed_by(self, subsumer: 'ConstrainedAtom') -> bool:
        """Return true if this c-atom (self) is independent of, or subsumed by, the subsumer.
        NOTE: This function is only as correct as '.subsumes' is."""
        return self.is_independent_from_other_clause(subsumer) or subsumer.subsumes(self)

    def subsumes(self, other: 'ConstrainedAtom') -> bool:
        """Does this atom (self) subsume the other atom (other)?
        We check this as they do in Forclift - checking if the atoms are split with respect to each other
        NOTE DEBUG: Trying simplified subsumes.
        TODO: work out if this is correct and if not, put back to original"""
        # quick check to handle the case where they're the same
        if self == other:
            return True

        mgu_eq_classes = self.get_constrained_atom_mgu_eq_classes(other)
        independent = mgu_eq_classes is None
        if independent:
            return False
        assert(mgu_eq_classes is not None)  # hack for type checking
        self_DNS_other = self.does_not_subsume(other, mgu_eq_classes)
        if not self_DNS_other:
            return True
        else: 
            return False

        # return self.needs_splitting(other) and not other.needs_splitting(self)

    def needs_splitting(self, other: 'ConstrainedAtom') -> bool:
        """Does this atom (self) need splitting with respect to the other atom (other)?
         Returns True if it does, and False otherwise."""
        mgu_eq_classes = self.get_constrained_atom_mgu_eq_classes(other)
        # check for independence
        if mgu_eq_classes is None:
            # print(f'independent')
            return False
        # check for subsumption
        else:
            # DEBUG
            # if not other.does_not_subsume(self, mgu_eq_classes):
            #     print(f'subsumes')
            return other.does_not_subsume(self, mgu_eq_classes)

    def does_not_subsume(self, other: 'ConstrainedAtom', mgu_eq_classes: 'EquivalenceClasses') -> bool:
        """Returns True if this ConstrainedAtom (self) does NOT subsume the
        other ConstrainedAtom (other) 
        NOTE: We apply the mgu as a substitution to the atoms so that they can
        be directly compared.  However, the EquivalenceClasses still refer to
        the pre-substitution atoms, so we mix and match.
        NOTE 2: I am going to treat free variables exactly like constants and see what happens.

        If any of the following are true, it does not subsume (the brackets are to avoid ambiguity):
        1) There is an mgu equivalence class (EC) that contains 
        (a constant or free variable) AND (a bound variable from self)
        2) There is an EC that contains TWO bound variables from self
        3) There is (an inequality between a bound variable and a constant) in the other that is
        not present in the self.
        4) There is (an inequality between bound variables in two ECs) in the other that is
        not present in the self.

        I am adding an experimental 5th rule to deal with domain variables
        (from getDomainShatteringMgu in Forclift)
        5) There is a bound variable in this atom (self) whose domain
        is a subset of the domain of a bound variable in the other atom (other)
        and the two variables are in the same EC.
        """

        # before everything else, rename to avoid overlapping bound variables
        other = other.make_variables_different(self)

        mgu_substitution = mgu_eq_classes.to_substitution()
        this_atom = self.substitute(mgu_substitution)
        other_atom = other.substitute(mgu_substitution)
        # DEBUG TODO: switch this back to returning True instead of numbers
        if any(( 
                len(eq_class.constants) > 0 
                or len(eq_class.variables.intersection(other.get_free_variables())) > 0
                or len(eq_class.variables.intersection(self.get_free_variables())) > 0
               )
               and len(eq_class.variables.intersection(other.bound_vars)) > 0
               for eq_class in mgu_eq_classes):
            # return "1"
            print("DNS 1")
            return True
        elif any(len(eq_class.variables.intersection(other.bound_vars)) >= 2
                 for eq_class in mgu_eq_classes):
            # return "2"
            print("DNS 2")
            return True
        elif any(inequality not in other_atom.get_constant_or_free_inequalities()
                 for inequality in this_atom.get_constant_or_free_inequalities()):
            # return "3"
            print("DNS 3")
            return True
        elif any(inequality not in other_atom.get_bound_variable_inequalities()
                 and inequality.is_not_trivial(this_atom)
                 for inequality in this_atom.get_bound_variable_inequalities()):
            # return "4"
            print("DNS 4")
            return True
        for eq_class in mgu_eq_classes:
            for term1 in eq_class:
                if term1 in self.bound_vars:
                    variable1 = term1
                    for term2 in eq_class:
                        if term2 in other.bound_vars:
                            variable2 = term2
                            this_domain = self.cs.get_domain_for(variable1)
                            other_domain = other.cs.get_domain_for(variable2)
                            if this_domain.is_strict_subset_of(other_domain):
                                print("DNS 5")
                                return True

        else:
            return False

    def make_variables_different(self, other_c_atom: 'ConstrainedAtom') -> 'ConstrainedAtom':
        """Make the bound variables of this c_atom (self) different to those in other_c_atom.
        NOTE: This is very similar to the method in UnitPropagation but serves a slightly different purpose"""

        overlapping_variables: List['LogicalVariable'] = []
        for variable in self.bound_vars:
            if variable in other_c_atom.bound_vars:
                overlapping_variables.append(variable)

        new_c_atom, new_other_c_atom = self, other_c_atom
        for variable in overlapping_variables:
            temp_cnf = CNF([new_c_atom, new_other_c_atom])  # taking advantage of existing methods in CNF
            sub_target = temp_cnf.get_new_logical_variable(variable.symbol)
            sub = Substitution([(variable, sub_target)])
            new_c_atom = new_c_atom.substitute(sub)
        return new_c_atom


