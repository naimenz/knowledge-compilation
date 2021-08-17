"""
Classes for literals, including the predicates and atoms that go into building literals.
"""

from kc.data_structures import Constant, LogicalVariable, EquivalenceClass, EquivalenceClasses
import typing

from typing import List, Sequence, Any, Set, Optional, FrozenSet, Iterable
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

    def is_smt(self) -> bool:
        """Is this an smt literal?"""
        # TODO: Figure out why mypy throws 'cannot determine type of predicate' here
        return isinstance(self.atom.predicate, SMTPredicate)  # type: ignore

    # TODO: refactor to sequence of if( and .. and ..) return True else False
    def __eq__(self, other: Any) -> bool: 
        """Two literals are equal if they have the same atom and polarity
        OR the predicate is SMT, one is positive and the other is negative, and the ranges work out"""
        if not isinstance(other, Literal):
            return False
        literals_identical = (self.atom == other.atom and self.polarity == other.polarity)
        if literals_identical:
            return True
        # handling the case of one-side unbounded smt atoms
        if self.is_smt() and other.is_smt():
            # TODO: Figure out why mypy throws 'cannot determine type of predicate' here
            self_pred, other_pred = self.atom.predicate, other.atom.predicate  # type: ignore
            if (self_pred.upper_bound == other_pred.lower_bound and self_pred.lower_bound == -other_pred.upper_bound == float('-inf')) \
                    or (self_pred.lower_bound == other_pred.upper_bound and self_pred.upper_bound == -other_pred.lower_bound == float('inf')):
                return True
        return False


        return isinstance(other, Literal) and self.atom == other.atom and self.polarity == other.polarity

    def __invert__(self) -> 'Literal':
        """Return a literal with the opposite polarity to this one.
        NOTE: We handle SMT literals differently if possible"""
        if self.is_smt():
            # TODO: Figure out why mypy can't determine type
            predicate = self.atom.predicate  # type: ignore
            if predicate.lower_bound == float('-inf'):
                new_predicate = SMTPredicate(predicate.name, predicate.arity, predicate.upper_bound, float('inf'))
                return Literal(Atom(new_predicate, self.atom.terms), True)  # type: ignore
            if predicate.upper_bound == float('inf'):
                new_predicate = SMTPredicate(predicate.name, predicate.arity, float('-inf'), predicate.lower_bound)
                return Literal(Atom(new_predicate, self.atom.terms), True)  # type: ignore
        return Literal(self.atom, not self.polarity)

    def __hash__(self) -> int:
        return hash((self.atom, self.polarity))

    def __str__(self) -> str:
        prefix = '¬' if self.polarity == False else ''
        return f'{prefix}{self.atom}'

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """The order is not important as long as it is consistent, so we will
        compare the order of the atom and the polarity"""
        if not isinstance(other, Literal): 
            raise NotImplementedError(f'Cannot compare Literal and {type(other)}')
        return (self.polarity, self.atom) < (other.polarity, other.atom)

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

    def is_smt(self) -> bool:
        """Is this an smt atom?"""
        return isinstance(self.predicate, SMTPredicate)

    def __eq__(self, other: Any) -> bool:
        """Two atoms are equal if they have the same predicate and the same terms"""
        return isinstance(other, Atom) and self.predicate == other.predicate and self.terms == other.terms

    def __hash__(self) -> int:
        return hash((self.predicate, self.terms))

    def __str__(self) -> str:
        term_strs = [str(term) for term in self.terms]
        return self.predicate.string_for_atom(term_strs)

    def string_for_wfomi(self) -> str:
        """For use in WFOMI, we need to give literals in a different format,
        e.g. 15 < age(X) < 99 as age_15_99(X)"""
        term_strs = [str(term) for term in self.terms]
        return self.predicate.string_for_wfomi(term_strs)

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """The order is not important as long as it is consistent, so we will
        compare the order of the predicate and the terms"""
        if not isinstance(other, Atom): 
            raise NotImplementedError(f'Cannot compare Atom and {type(other)}')
        return (self.predicate, self.terms) < (other.predicate, other.terms)



# class GroundAtom(Atom):
#     """
#     A class for ground FOL atoms.
#     These are the same as regular atoms but all terms are constants
#     """
#     def __init__(self, predicate: 'Predicate', terms: List['Constant']) -> None:
#         cast_terms = cast(List['LogicalTerm'], terms) # hack for type checking
#         super().__init__(predicate, cast_terms)

#     @staticmethod
#     def build_from_atom_substitution(atom: 'Atom', substitution: 'Substitution') -> 'GroundAtom':
#         """Construct a ground atom from a substitution to a non-ground atom

#         NOTE: assumes that all variables in the atom appear in the substitution, or it'll throw a key error.
#         TODO: should this function be here, or in Atom, or stand-alone?"""
#         ground_terms: List['Constant'] = []
#         for term in atom.terms:
#             if isinstance(term, Constant):
#                 ground_terms.append(term)
#             else:
#                 variable = cast(LogicalVariable, term) # hack for type checking 
#                 sub = substitution[variable]
#                 if not isinstance(sub, Constant):
#                     raise ValueError('Substitution for ground atom should be closing (i.e. no variables)')
#                 ground_terms.append(sub)
#         return GroundAtom(atom.predicate, ground_terms)


class Predicate:
    """
    A class for FOL predicates.
    This consists of a predicate name and an arity (i.e. how many arguments)
    """

    def __init__(self, name: str, arity: int) -> None:
        self.name = name
        self.arity = arity

    def is_smt(self) -> bool:
        """Is this an smt predicate?"""
        return isinstance(self, SMTPredicate)

    def __eq__(self, other: Any) -> bool:
        """Two predicates are equal if they have the same name and the same arity
        NOTE: We do not need to handle the SMTPredicate case, since the way python handles
        __eq__, the method of the other instance will be called even if it is on the right"""
        return isinstance(other, Predicate) and self.name == other.name and self.arity == other.arity

    def __hash__(self) -> int:
        return hash((self.name, self.arity))

    def __str__(self) -> str:
        return f"{self.name}/{self.arity}"

    def string_for_atom(self, term_strs: Iterable[str]) -> str:
        """For representation inside an atom, we need to provide a different function,
        with space to include terms"""
        return f"{self.name}({','.join(term_strs)})"

    def __repr__(self) -> str:
        return self.__str__()

    def __lt__(self, other: Any) -> bool:
        """The order is not important as long as it is consistent, so we will
        compare the order of the symbol and the arity"""
        if not isinstance(other, Predicate): 
            raise NotImplementedError(f'Cannot compare Predicate and {type(other)}')
        return (self.name, self.arity) < (other.name, other.arity)

class SMTPredicate(Predicate):
    """
    A class for SMT predicates.
    This consists of a predicate name (which is really a function symbol name for SMT predicates),
    an arity, AND upper and lower bounds (either or both of which can be inf/-inf)
    """
    def __init__(self, name: str, arity: int, lower_bound: float, upper_bound: float) -> None:
        super(SMTPredicate, self).__init__(name, arity)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def is_smt(self) -> bool:
        """Is this an smt predicate?"""
        return True

    def __eq__(self, other: Any) -> bool:
        """Two predicates are equal if they have the same name and the same arity"""
        return isinstance(other, SMTPredicate) and self.name == other.name and self.arity == other.arity and \
                self.lower_bound == other.lower_bound and self.upper_bound == other.upper_bound
    
    def __hash__(self) -> int:
        return hash((self.name, self.arity, self.lower_bound, self.upper_bound))

    def __str__(self) -> str:
        le_string = ' \u2264 '
        return f"{self.lower_bound}{le_string}{self.name}/{self.arity} < {self.upper_bound}"

    def string_for_atom(self, term_strs: Iterable[str]) -> str:
        """For representation inside an atom, we need to provide a different function,
        with space to include terms"""
        le_string = ' \u2264 '
        return f"[{self.lower_bound}{le_string}{self.name}({','.join(term_strs)}) < {self.upper_bound}]"

    def string_for_wfomi(self, term_strs: Iterable[str]) -> str:
        """For use in WFOMI, we need to give literals in a different format,
        e.g. 15 < age(X) < 99 as age_15_99(X)
        TODO: Decide on how many decimal places"""
        return f"{self.name}_{self.lower_bound:.0f}_{self.upper_bound:.0f}({','.join(term_strs)})"

    def __lt__(self, other: Any) -> bool:
        """The order is not important as long as it is consistent, so we will
        compare the order of the symbol and the arity"""
        if not isinstance(other, Predicate): 
            raise NotImplementedError(f'Cannot compare Predicate and {type(other)}')
        return (self.name, self.arity, self.lower_bound, self.upper_bound) < (other.name, other.arity, self.lower_bound, self.upper_bound)
