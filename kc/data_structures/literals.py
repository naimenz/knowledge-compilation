"""
Classes for literals, including the predicates and atoms that go into building literals.
"""

from kc.data_structures import Constant, LogicalVariable, EquivalenceClass, EquivalenceClasses
import typing

from typing import List, Sequence, Any, Set, Optional, FrozenSet
from typing import cast
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import Substitution, LogicalTerm


class Literal:
    """
    A class for FOL literals.
    This consists of an atom with a polarity (i.e. is it negated?)
    """

    def __init__(self, atom: 'Atom', polarity: bool = True) -> None:
        self.atom = atom
        self.polarity = polarity

    def substitute(self, substitution: 'Substitution') -> 'Literal':
        """Return a new Literal, the result of applying substitution to the current Literal"""
        return Literal(self.atom.substitute(substitution), self.polarity)

    @property
    def variables(self) -> FrozenSet['LogicalVariable']:
        """Return only the variables in the terms of this literal"""
        return self.atom.variables

    @property
    def constants(self) -> FrozenSet['Constant']:
        """Return only the constants in the terms of this literal"""
        return self.atom.constants

    # TODO: refactor to sequence of if( and .. and ..) return True else False
    def __eq__(self, other: Any) -> bool: 
        """Two literals are equal if they have the same atom and polarity"""
        return isinstance(other, Literal) and self.atom == other.atom and self.polarity == other.polarity

    def __invert__(self) -> 'Literal':
        """Return a literal with the opposite polarity to this one."""
        return Literal(self.atom, not self.polarity)

    def __hash__(self) -> int:
        return hash((self.atom, self.polarity))

    def __str__(self) -> str:
        prefix = '¬' if self.polarity == False else ''
        return f'{prefix}{self.atom}'

    def __repr__(self) -> str:
        return self.__str__()

class Atom:
    """
    A class for FOL atoms.
    This consists of a predicate and a tuple of terms.
    """

    def __init__(self, predicate: 'Predicate', terms: Sequence['LogicalTerm']) -> None:
        assert(len(terms) == predicate.arity)
        self.predicate = predicate
        self.terms = tuple(terms)
        variables, constants = set(), set()
        for term in terms:
            if isinstance(term, LogicalVariable):
                variables.add(term)
            elif isinstance(term, Constant):
                constants.add(term)
        self._variables = frozenset(variables)
        self._constants = frozenset(constants)

    def substitute(self, substitution: 'Substitution') -> 'Atom':
        """Return a new Atom, the result of applying substitution to the current Atom."""
        new_terms: List['LogicalTerm'] = []
        for term in self.terms:
            sub_value = substitution[term]
            if not sub_value is None:
                new_term = sub_value
            else:
                new_term = term
            new_terms.append(new_term)
        return Atom(self.predicate, new_terms)

    @property
    def variables(self) -> FrozenSet['LogicalVariable']:
        """Return only the logical variables that appear in the terms of this atom"""
        return self._variables


    @property
    def constants(self) -> FrozenSet['Constant']:
        """Return only the constants that appear in the terms of this atom"""
        return self._constants

    def get_unconstrained_atom_mgu_eq_classes(self: 'Atom', other_atom: 'Atom') -> Optional['EquivalenceClasses']:
        """Compute the mgu of this atom with another using equivalence classes,

        Returns a list of equivalence classes or None, if unsuccessful.
        """
        if self.predicate != other_atom.predicate:
            return None

        # we build up equivalence classes, starting from equivalences between terms in the same position
        # of the two atoms
        term_pairs = zip(self.terms, other_atom.terms)
        # NOTE: Including singleton equivalence classes so they aren't missed by ISG
        # # we only include equivalence classes with more than one unique element
        initial_eq_classes = EquivalenceClasses( EquivalenceClass([t1, t2]) for t1, t2 in term_pairs)# if t1 != t2 )

        final_eq_classes = initial_eq_classes.make_self_consistent()
        if any(eq_class.is_inconsistent() for eq_class in final_eq_classes):
            return None
        else:
            return final_eq_classes

    def get_unconstrained_atom_mgu_substitution(self: 'Atom', other_atom: 'Atom') -> Optional['Substitution']:
        """Compute the mgu of this atom with another using equivalence classes,
        as done in Forclift. Before returning, I convert the equivalence classes into a
        substitution.

        Returns None if the mgu doesn't exist (the atoms are independent), and a Substitution if it does.
        """
        eq_classes = self.get_unconstrained_atom_mgu_eq_classes(other_atom)
        if not eq_classes is None:
            return eq_classes.to_substitution()
        else:
            return None

    def __eq__(self, other: Any) -> bool:
        """Two atoms are equal if they have the same predicate and the same terms"""
        return isinstance(other, Atom) and self.predicate == other.predicate and self.terms == other.terms

    def __hash__(self) -> int:
        return hash((self.predicate, self.terms))

    def __str__(self) -> str:
        term_strs = [str(term) for term in self.terms]
        return f"{self.predicate.name}({', '.join(term_strs)})"

    def __repr__(self) -> str:
        return self.__str__()


class GroundAtom(Atom):
    """
    A class for ground FOL atoms.
    These are the same as regular atoms but all terms are constants
    """
    def __init__(self, predicate: 'Predicate', terms: List['Constant']) -> None:
        cast_terms = cast(List['LogicalTerm'], terms) # hack for type checking
        super().__init__(predicate, cast_terms)

    @staticmethod
    def build_from_atom_substitution(atom: 'Atom', substitution: 'Substitution') -> 'GroundAtom':
        """Construct a ground atom from a substitution to a non-ground atom

        NOTE: assumes that all variables in the atom appear in the substitution, or it'll throw a key error.
        TODO: should this function be here, or in Atom, or stand-alone?"""
        ground_terms: List['Constant'] = []
        for term in atom.terms:
            if isinstance(term, Constant):
                ground_terms.append(term)
            else:
                variable = cast(LogicalVariable, term) # hack for type checking 
                sub = substitution[variable]
                if not isinstance(sub, Constant):
                    raise ValueError('Substitution for ground atom should be closing (i.e. no variables)')
                ground_terms.append(sub)
        return GroundAtom(atom.predicate, ground_terms)


class Predicate:
    """
    A class for FOL predicates.
    This consists of a predicate name and an arity (i.e. how many arguments)
    """

    def __init__(self, name: str, arity: int) -> None:
        self.name = name
        self.arity = arity

    def __eq__(self, other: Any) -> bool:
        """Two predicates are equal if they have the same name and the same arity"""
        return isinstance(other, Predicate) and self.name == other.name and self.arity == other.arity

    def __hash__(self) -> int:
        return hash((self.name, self.arity))

    def __str__(self) -> str:
        return f"{self.name}/{self.arity}"

    def __repr__(self) -> str:
        return self.__str__()

