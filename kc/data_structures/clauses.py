"""
Classes for clauses in FOL-DC.
This includes constrained AND unconstrained clauses.

TODO: Figure out if the inheritance structure for UnitClause and ConstrainedAtom makes sense.
"""

from kc.data_structures import Literal, Atom, LogicalVariable, Constant, EquivalenceClass
from kc.data_structures import ConstraintSet, EquivalenceClasses

from functools import reduce
from abc import ABC, abstractmethod

from typing import List, TypeVar, Iterable, Any, Sequence, Set, Optional
from typing import cast


class Clause(ABC):
    """Abstract base class for constrained and unconstrained clauses"""
    # @abstractmethod
    # def apply_substitution(self: 'Clause', substitution: 'Substitution') -> 'Clause':
    #     pass


class UnconstrainedClause(Clause):
    """An FOL unconstrained clause.
    This consists of a set of FOL literals, which form a disjunction
    """

    def __init__(self, literals: Iterable['Literal']) -> None:
        """NOTE: using a set to represent literals.
        This is justified because the order is unimportant and repeated literals 
        in a clause are redundant (they cannot change the disjunction)"""
        self.literals = frozenset(literals)

    # def apply_substitution(self: 'UnconstrainedClause', substitution: 'Substitution') -> 'UnconstrainedClause':
    #     """Return a new UnconstrainedClause, the result of applying substitution to this UnconstrainedClause"""
    #     return self.__class__(literal.apply_substitution(substitution) for literal in self.literals)

    def __eq__(self, other: Any) -> bool:
        """Two unconstrained clauses are equal if they have the same literals

        NOTE: for now the ordering has to be the same
        TODO: make literals hashable so I can compare as sets"""
        if not isinstance(other, UnconstrainedClause):
            return False
        if len(self.literals) != len(other.literals): 
            return False
        same_literals = all(self_l == other_l for self_l, other_l in zip(self.literals, other.literals))
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
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        self.unconstrained_clause = unconstrained_clause
        self.bound_vars = frozenset(bound_vars)
        self.cs = cs

    # def apply_substitution(self: 'C', substitution: 'Substitution') -> 'C':
    #     """Return a new ConstrainedClause, the result of applying substitution to this ConstrainedClause
    #     NOTE: assumes that the bound vars aren't substituted"""
    #     new_unconstrained_clause = self.unconstrained_clause.apply_substitution(substitution)
    #     new_cs = self.cs.apply_substitution(substitution)
    #     return self.__class__(new_unconstrained_clause, self.bound_vars, new_cs)

    @property
    def all_literal_variables(self) -> Set['LogicalVariable']:
        """Return all variables that appear in any literal of the clause"""
        all_vars: Set['LogicalVariable'] = set()
        for literal in self.unconstrained_clause.literals:
            all_vars = all_vars.union(literal.variables)
        return all_vars

    @property
    def all_variables(self) -> Set['LogicalVariable']:
        """Return ALL variables, regardless of whether they are bound, in the cs,
        or in the literals"""
        literal_vars = self.all_literal_variables
        bound_vars = self.bound_vars
        constraint_vars = self.constraint_variables 
        return literal_vars.union(bound_vars, constraint_vars)

    @property
    def root_variables(self) -> Set['LogicalVariable']:
        """Return the root variables of this clause 
        (i.e. the variables that appear in every literal of the clause)"""
        # first, collect ALL variables that appear in any literal
        root_vars = self.all_literal_variables
        # then get only the ones that appear in every literal
        for literal in self.unconstrained_clause.literals:
            root_vars = root_vars.intersection(literal.variables)
        return root_vars

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

    def clauses_independent(self, other_clause: 'ConstrainedClause') -> bool:
        """Is this clause independent of the other clause?"""
        if self._is_conditional_contradiction() or other_clause._is_conditional_contradiction():
            return True

        all_independent = True
        for c_atom in self.get_constrained_atoms():
            for other_c_atom in other_clause.get_constrained_atoms():
                if c_atom.constrained_atoms_unify(other_c_atom):
                    all_independent = False
        return all_independent

    def _is_conditional_contradiction(self) -> bool:
        """This name is from Forclift.
        I think this means that the clause contains no literals and so
        its grounding is empty, meaning it is independent of everything."""
        return len(self.unconstrained_clause.literals) == 0

    def get_constrained_atoms(self) -> List['ConstrainedAtom']:
        """For this clause, return a list of all the constrained atoms in the clause"""
        constrained_atoms = []
        for literal in self.unconstrained_clause.literals:
            constrained_atom = self._build_constrained_atom_from_literal(literal)
            constrained_atoms.append(constrained_atom)
        return constrained_atoms

    def _build_constrained_atom_from_literal(self, literal: 'Literal') -> 'ConstrainedAtom':
        """Build a constrained atom given a specific literal using the bound variables
        and constraints of this clause"""
        atom = literal.atom 
        positive_literal = Literal(atom, True)
        # constrained atoms subclass ConstrainedClause, so need to be built of an
        # unconstrained clause
        unconstrained_atom = UnconstrainedClause([positive_literal])
        constrained_atom = ConstrainedAtom(unconstrained_atom, self.bound_vars, self.cs)
        return constrained_atom

    def __eq__(self, other: Any) -> bool:
        """Two constrained literals are equal if they have the same unconstrained literals, the same constraint sets,
         and the same bound variables"""
        if not isinstance(other, ConstrainedClause):
            return False
        same_u_clause = (self.unconstrained_clause == other.unconstrained_clause)
        if len(self.bound_vars) != len(other.bound_vars): 
            return False
        same_bound_vars = all(self_v == other_v for self_v, other_v in zip(self.bound_vars, other.bound_vars))
        same_cs = (self.cs == other.cs)
        return same_u_clause and same_bound_vars and same_cs

    def __hash__(self) -> int:
       return hash((self.unconstrained_clause, self.bound_vars, self.cs))

    def __str__(self) -> str:
        bound_vars_strs = [str(var) for var in self.bound_vars]
        for_all_string = '\u2200'
        return f"{for_all_string}{{{', '.join(bound_vars_strs)}}}, {self.cs} : {self.unconstrained_clause}"

    def __repr__(self) -> str:
        return self.__str__()


class UnitClause(ConstrainedClause):
    """AN FOL-DC unit clause.
    This is a constrained clause with a single literal
    """
    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(unconstrained_clause.literals) == 1) # ensure that this is a unit clause
        self.literal = list(unconstrained_clause.literals)[0] # convert to 1-item list and get item
        super(UnitClause, self).__init__(unconstrained_clause, bound_vars, cs)


class ConstrainedAtom(UnitClause):
    """An FOL-DC constrained atom.
    This is a constrained clause with a single atom (positive literal).
    """
    def __init__(self,
            unconstrained_clause: 'UnconstrainedClause',
            bound_vars: Iterable['LogicalVariable'],
            cs: 'ConstraintSet') -> None:
        assert(len(unconstrained_clause.literals) == 1) # ensure that this is a unit clause
        super(ConstrainedAtom, self).__init__(unconstrained_clause, bound_vars, cs)
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


# Type variable for arbitrary clauses so I can reuse apply_substitution
# ideally this would be defined above but I think mypy is bugged
C = TypeVar('C', bound='ConstrainedClause') 

