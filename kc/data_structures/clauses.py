"""
Classes for clauses in FOL-DC.
This includes constrained AND unconstrained clauses.

TODO: Figure out if the inheritance structure for UnitClause and ConstrainedAtom makes sense.
"""

from kc.data_structures import Literal, Atom, LogicalVariable, Constant, ConstraintSet

from functools import reduce
from abc import ABC, abstractmethod

from typing import List, TypeVar, Iterable, Any, Sequence, Set, Optional
from typing import cast 
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import Atom, Substitution, ConstraintSet, EquivalenceClasses

# Type variable for arbitrary clauses so I can reuse apply_substitution
C = TypeVar('C', bound='ConstrainedClause') 

class Clause(ABC):
    """Abstract base class for constrained and unconstrained clauses"""
    def __init__(self, literals: Sequence['Literal']) -> None:
        self.literals = frozenset(literals)

    @abstractmethod
    def apply_substitution(self: 'Clause', substitution: 'Substitution') -> 'Clause':
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

    def clauses_independent(self, other_clause: 'Clause') -> bool:
        """Is this clause independent of the other clause?"""
        if self.has_no_literals() or other_clause.has_no_literals():
            return True

        all_independent = True
        for c_atom in self.get_constrained_atoms():
            for other_c_atom in other_clause.get_constrained_atoms():
                if c_atom.constrained_atoms_unify(other_c_atom):
                    all_independent = False
        return all_independent

    # TODO: this may be wrong, may need every literal to subsume another
    def is_subsumed_by_clause(self, subsumer: 'Clause') -> bool:
        """Returns true if this Clause is subsumed by the other.
        We do this by looking at constrained atoms -- if ANY literal in the subsumer
        subsumes ANY literal in self, then it is subsumed (because disjunction).
        NOTE: This is only as correct as 'is_subsumed_by_c_atom' is."""
        if self.has_no_literals():
            return True

        for c_literal in self.get_constrained_literals():
            for other_c_literal in subsumer.get_constrained_literals():
                if c_literal.is_subsumed_by_literal(other_c_literal):
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

    def __init__(self, literals: Iterable['Literal']) -> None:
        """NOTE: using a set to represent literals.
        This is justified because the order is unimportant and repeated literals 
        in a clause are redundant (they cannot change the disjunction)"""
        self.literals = frozenset(literals)

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

    # TODO:  maybe change name to substitute?
    def apply_substitution(self: 'UnconstrainedClause', substitution: 'Substitution') -> 'UnconstrainedClause':
        """Return a new UnconstrainedClause, the result of applying substitution to this UnconstrainedClause"""
        return self.__class__(literal.apply_substitution(substitution) for literal in self.literals)

    def __eq__(self, other: Any) -> bool:
        """Two unconstrained clauses are equal if they have the same literals

        NOTE: for now the ordering has to be the same
        TODO: make literals hashable so I can compare as sets"""
        if not isinstance(other, UnconstrainedClause):
            return False
        if len(self.literals) != len(other.literals): 
            return False
        same_literals = (self.literals == other.literals)
        return same_literals

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

    def apply_substitution(self: 'C', substitution: 'Substitution') -> 'C':
        """Return a new ConstrainedClause, the result of applying substitution to this ConstrainedClause
        NOTE: assumes that the bound vars aren't substituted"""
        new_literals = [literal.apply_substitution(substitution) for literal in self.literals]
        new_cs = self.cs.apply_substitution(substitution)
        new_bound_vars = [substitution[var] for var in self.bound_vars]
        assert(all(isinstance(term, LogicalVariable) for term in new_bound_vars)) # should not have constants in bound vars
        bound_vars = cast(List['LogicalVariable'], new_bound_vars) # hack for type checking
        return self.__class__(new_literals, bound_vars, new_cs)


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


    def clauses_independent(self, other_clause: 'Clause') -> bool:
        """Is this clause independent of the other clause?"""
        if self.has_no_literals() or other_clause.has_no_literals():
            return True

        all_independent = True
        for c_atom in self.get_constrained_atoms():
            for other_c_atom in other_clause.get_constrained_atoms():
                if c_atom.constrained_atoms_unify(other_c_atom):
                    all_independent = False
        return all_independent

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
        if not isinstance(other, ConstrainedClause):
            return False
        same_literals = (self.literals == other.literals)
        if len(self.bound_vars) != len(other.bound_vars): 
            return False
        same_bound_vars = (self.bound_vars == other.bound_vars)
        same_cs = (self.cs == other.cs)
        return same_literals and same_bound_vars and same_cs

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
          and self.to_c_atom().is_subsumed_by_c_atom(subsumer.to_c_atom()):
                return True
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
        NOTE: This function is only as correct as 'is_subsumed_by_c_atom'."""
        return self.clauses_independent(subsumer) or self.is_subsumed_by_c_atom(subsumer)

    # TODO: Continue improving this function to work in more cases.
    def is_subsumed_by_c_atom(self, subsumer: 'ConstrainedAtom') -> bool:
        """Does the subsumer (A) subsume this c_atom, the subsumed (B)?
        NOTE: This is a work in progress, and currently uses the following (incomplete) rules:
        1) A and B must unify, producing equivalence classes.
        2) Each equivalence class between a variable X in A and a constant c in B must have
        c in the shared domain for X (after processing inequality constraints)
        3) For each equivalence class, its shared domain in B must be a subset of its shared domain in A.
        4) Each equivalence class must contain only ONE term from B.
        5) FOR NOW: A free variable anywhere in A breaks subsumption. 
        A free variable in the constraint set of B does not, but a free variable in its atom does. 
        (THIS IS WRONG, BUT IS A FIRST DRAFT)
        The main change I want to make is to check if both clauses contain equivalent free
        variables, but checking this is hard, especially for the constraint set
        """
        # aliases because this is how I've been using them in my notes
        A, B = subsumer, self
        # TODO: remove this hack of checking for equality, which only works some of the time
        # if A == B:
        #     return True
        eq_classes = A.get_constrained_atom_mgu_eq_classes(B)
        # 1)
        if eq_classes is None:
            print("DEBUG: Didn't unify")
            return False
        for eq_class in eq_classes:
            var_eq_class = eq_class.get_variables_only()
            A_shared_domain = var_eq_class.get_shared_domain_from_cs(A.cs)
            B_shared_domain = var_eq_class.get_shared_domain_from_cs(B.cs)

            # 2) - this is a long, potentially slow check TODO: make it neater
            for eq_term in eq_class:
                if isinstance(eq_term, Constant):
                    for i, arg_term in enumerate(B.atom.terms):
                        if eq_term == arg_term and isinstance(A.atom.terms[i], LogicalVariable):
                            if not eq_term in A_shared_domain:
                                print(f'DEBUG: Constant {eq_term} not in {A_shared_domain=}')
                                return False

            # 3)
            if not B_shared_domain.is_subset_of(A_shared_domain):
                print(f'DEBUG: {B_shared_domain=} is not a subset of {A_shared_domain=}')
                return False

            # 4)
            if len(eq_class.members.intersection(B.literal_variables)) > 1:
                print(f'DEBUG: {eq_class=} and {B.literal_variables=} overlap in more than one place')
                return False

        # 5)
        A_free_variables = A.all_variables.difference(A.bound_vars)
        if len(A_free_variables) > 0:
            print(f'DEBUG: {A_free_variables=}')
            return False

        B_literal_free_variables = B.literal_variables.difference(B.bound_vars)
        if len(B_literal_free_variables) > 0:
            print(f'DEBUG: {B_literal_free_variables=}')
            print(B.literal_variables)
            print(B.bound_vars)
            return False

        return True



