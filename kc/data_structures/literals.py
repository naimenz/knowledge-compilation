"""
Classes for literals, including the predicates and atoms that go into building literals.
"""

from kc.data_structures.logicalterms import *
from kc.data_structures.substitutions import *
import typing

from typing import List, Sequence, Any 
from typing import cast

class Literal:
    """
    A class for FOL literals.
    This consists of an atom with a polarity (i.e. is it negated?)
    """

    def __init__(self, atom: 'Atom', polarity: bool) -> None:
        self.atom = atom
        self.polarity = polarity

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

    def __init__(self, predicate: 'Predicate', terms: List['LogicalTerm']) -> None:
        assert(len(terms) == predicate.arity)
        self.predicate = predicate
        self.terms = terms

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
                ground_terms.append(substitution[variable])
        return GroundAtom(atom.predicate, ground_terms)

    def __str__(self) -> str:
        term_strs = [str(term) for term in self.terms]
        return f"{self.predicate.name}({', '.join(term_strs)})"

    def __repr__(self) -> str:
        return self.__str__()


class Predicate:
    """
    A class for FOL predicates.
    This consists of a predicate name and an arity (i.e. how many arguments)
    """

    def __init__(self, name: str, arity: int) -> None:
        self.name = name
        self.arity = arity

    def __str__(self) -> str:
        return f"{self.name}/{self.arity}"

    def __repr__(self) -> str:
        return self.__str__()

if __name__ == '__main__':
    pred = Predicate('smokes', 4)
    print(pred)

    terms: List[Any] = [Constant('a'), LogicalVariable('X'), Constant('b'), LogicalVariable('Y')]
    atom = Atom(pred, terms)
    print(atom)

    ground_atom = Atom(pred, [terms[0], terms[2], terms[0], terms[0]])
    print("ground atom",ground_atom)

    substitution = Substitution([(terms[1], terms[0]), (terms[3], terms[2])])
    built_atom = GroundAtom.build_from_atom_substitution(atom, substitution)
    print("built atom", built_atom)

    literal = Literal(atom, False)
    print(literal)



