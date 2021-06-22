
"""
Classes for types of domain terms, which are sets of constants and domain variables
"""

from abc import ABC, abstractmethod
from kc.data_structures.logicalterms import *

from typing import List, Set, Any, FrozenSet, Iterable
from typing import cast, TypeVar

# generic type for domain term or its subclasses
TDomain = TypeVar('TDomain', bound='DomainTerm')

class DomainTerm(ABC):
    """
    An abstract base class for domain terms.
    Terms are either sets of constants or domain variables, so can never be instantiated directly.
    """

    @staticmethod
    def union(*args: 'DomainTerm') -> 'DomainTerm':
        """Experimental function for taking the union of domain terms.
        NOTE: at the moment we only deal with SetOfConstants"""
        # TODO: maybe return a special EmptyDomain?
        if len(args) == 0:
            return SetOfConstants([])
        if all(isinstance(term, SetOfConstants) for term in args):
            domains = cast(List['SetOfConstants'], args) # hack for type checking
            constants: Set['Constant'] = set.union(*[set(domain.constants) for domain in domains])
            return SetOfConstants(constants)

        else:
            raise NotImplementedError('union only works for SetOfConstants for now')

    @staticmethod
    def intersection(*args: 'DomainTerm') -> 'DomainTerm':
        """Experimental function for taking the intersection of domain terms.
        NOTE: at the moment we only deal with SetOfConstants"""
        # TODO: maybe return a special EmptyDomain?
        if len(args) == 0:
            return SetOfConstants([])
        if all(isinstance(term, SetOfConstants) for term in args):
            domains = cast(List['SetOfConstants'], args) # hack for type checking
            constants: Set['Constant'] = set.intersection(*[set(domain.constants) for domain in domains])
            return SetOfConstants(constants)
        else:
            raise NotImplementedError('intersection only works for SetOfConstants for now')

    def is_subset_of(self, other: 'DomainTerm') -> bool:
        """Return True if THIS SetOfConstants (self) is a subset of 'other' and False otherwise
        NOTE: For now only works with SetOfConstants"""
        if all(isinstance(term, SetOfConstants) for term in (self, other) ):
            domains = cast(List['SetOfConstants'], (self, other) ) # hack for type checking
            constants: Set['Constant'] = set.intersection(*[set(domain.constants) for domain in domains])
            return SetOfConstants(constants) == self
        else:
            raise NotImplementedError('is_subset_of only works for SetOfConstants for now')

    def is_superset_of(self, other: 'DomainTerm') -> bool:
        """Return True if THIS SetOfConstants (self) is a superset of 'other' and False otherwise
        NOTE: For now only works with SetOfConstants"""
        if all(isinstance(term, SetOfConstants) for term in (self, other) ):
            domains = cast(List['SetOfConstants'], (self, other) ) # hack for type checking
            constants: Set['Constant'] = set.intersection(*[set(domain.constants) for domain in domains])
            return SetOfConstants(constants) == other
        else:
            raise NotImplementedError('is_superset_of only works for SetOfConstants for now')
        

    @property
    @abstractmethod
    def size(self) -> int:
        """I think domain variables will need sizes too"""


class SetOfConstants(DomainTerm):
    """
    A set of FOL constants. 
    """
    def __init__(self, constants: Iterable['Constant']) -> None:
        self._constants = frozenset(constants)

    @property
    def constants(self) -> FrozenSet['Constant']:
        """We use a property because the sets of constants should
        never be changed once set."""
        return self._constants

    def difference(self, other: 'SetOfConstants') -> 'SetOfConstants':
        """Return the set difference between this set of constants and the other
        as a SetOfConstants"""
        new_constants = self.constants - other.constants
        return SetOfConstants(new_constants)

    @property
    def size(self) -> int:
        return len(self.constants)

    def __iter__(self) -> Iterable['Constant']:
        """Replacing __contains__ with iter for more flexibility"""
        return self.constants

    def __contains__(self, other: Any) -> bool:
        """Returns True if this SetOfConstants contains the other thing."""
        if not isinstance(other, 'Constant'):
            raise ValueError(f'Only constants can be in SetOfConstants, not {type(other)}')
        return other in self.constants


    def __eq__(self, other: Any) -> bool:
        """Two sets of constants are equal if their constants are the same"""
        if not isinstance(other, SetOfConstants):
            return False
        return self.constants == other.constants

    def __hash__(self) -> int:
        return hash(self.constants)

    def __str__(self) -> str:
        constant_strs = [str(constant) for constant in self.constants]
        return f"{{{', '.join(constant_strs)}}}"

    def __repr__(self) -> str:
        return self.__str__()

class DomainVariable(DomainTerm):
    """
    A FOL domain variable.
    """
    def __init__(self, symbol: str) -> None:
        self.symbol = symbol

    def __eq__(self, other: Any) -> bool:
        """Two domain variables are equal if they have the same symbol.

        NOTE: I'm not sure if this is the right call, since symbols are sometimes reused.
        TODO: figure out what to do about symbol reuse and remove warning."""
        if not isinstance(other, DomainVariable):
            return False
        # this is the case I'm not sure about
        if (self.symbol == other.symbol) and not (self is other):
            print(f"WARNING: Treating DomainVariables {self.symbol} as equal but are not the same object")
            return True
        elif (self is other):
            return True
        else:
            return False

    def __hash__(self) -> int:
        return hash(self.symbol)

    def __str__(self) -> str:
        return f'{self.symbol}'

    def __repr__(self) -> str:
        return self.__str__()

