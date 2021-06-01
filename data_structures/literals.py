"""
Classes for literals, including the predicates and atoms that go into building literals.
"""

from logicalterms import *
import typing

class Literal:
    """
    A class for FOL literals.
    This consists of an atom with a polarity (i.e. is it negated?)
    """

    def __init__(self, atom: 'Atom', polarity: bool) -> None:
        self.atom = atom
        self.polarity = polarity

    def __str__(self) -> str:
        prefix = '¬' if not self.polarity else ''
        return f'{prefix}{self.atom}'


class Atom:
    """
    A class for FOL atoms.
    This consists of a predicate and a tuple of terms.
    """

    def __init__(self, predicate: 'Predicate', terms: list['LogicalTerm']) -> None:
        assert(len(terms) == predicate.arity)
        self.predicate = predicate
        self.terms = terms

    def __str__(self) -> str:
        term_strs = [str(term) for term in self.terms]
        return f"{self.predicate.name}({', '.join(term_strs)})"


class Predicate:
    """
    A class for FOL predicates.
    This consists of a predicate name and an arity (i.e. how many arguments)
    """

    def __init__(self, name: str, arity: int) -> None:
        self.name = name
        self.arity = arity

    def __str__(self) -> str:
        return f"{self.name}({self.arity})"


if __name__ == '__main__':
    pred = Predicate('smokes', 4)
    print(pred)

    terms = [Constant('a'), LogicalVariable('X'), Constant('b'), LogicalVariable('Y')]
    atom = Atom(pred, terms)
    print(atom)

    literal = Literal(atom, False)
    print(literal)

