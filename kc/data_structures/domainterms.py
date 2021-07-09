
"""
Classes for types of domain terms, which are sets of constants and domain variables
"""

from abc import ABC, abstractmethod
from kc.data_structures import Constant

from typing import List, Set, Any, FrozenSet, Iterable
from typing import cast

class DomainTerm(ABC):
    """
    An abstract base class for domain terms.
    Terms are either sets of constants or domain variables, so can never be instantiated directly.
    TODO: Improve DomainVariable
    """

    @abstractmethod
    def difference(self, other: 'DomainTerm') -> 'DomainTerm':
        """Get the difference between this domain term and another"""
        pass

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

    def difference(self, other: 'DomainTerm') -> 'DomainTerm':
        """Return the set difference between this set of constants and the other
        as a SetOfConstants (if both are SetOfConstants)"""
        if not isinstance(other, SetOfConstants):
            raise NotImplementedError('Cannot compute difference of SetOfConstants with DomainVariable')
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
        if not isinstance(other, Constant):
            raise ValueError(f'Only constants can be in SetOfConstants, not {type(other)}')
        return other in self.constants


    def __eq__(self, other: Any) -> bool:
        """Two sets of constants are equal if their constants are the same"""
        return isinstance(other, SetOfConstants) and self.constants == other.constants

    def __hash__(self) -> int:
        return hash(self.constants)

    def __str__(self) -> str:
        constant_strs = [str(constant) for constant in self.constants]
        return f"{{{', '.join(constant_strs)}}}"

    def __repr__(self) -> str:
        return self.__str__()

class ProperDomain:
    """This is supposed to be an abstract class representing domain terms
    that can be the proper domain of a logical variable rather than just some constants.
    TODO: Work out if this approach makes sense"""

class RootDomain(SetOfConstants, ProperDomain):
    """This is a class to represent a RootDomain, i.e. a specific set of constants
    given as part of the user-defined input. New RootDomains will not be created
    during compilation."""
    def __init__(self, constants: Iterable['Constant'], symbol: str) -> None:
        super(RootDomain, self).__init__(constants)
        self.symbol = symbol
        self.children: List['DomainTerm'] = []  # will be updated as children are created

    def __eq__(self, other: Any) -> bool:
        """Since RootDomains should not be created except at the start,
        there should only ever one of each. We can just check if the symbols and constants match."""
        return isinstance(other, RootDomain) and self.constants == other.constants and self.symbol == other.symbol

    def __hash__(self) -> int:
        return hash((self.symbol, self.constants))

    def __str__(self) -> str:
        return f"{self.symbol}"

    def __repr__(self) -> str:
        return self.__str__()

class DomainVariable(DomainTerm, ProperDomain):
    """
    A FOL domain variable.
    NOTE: This is based quite heavily on Forclift's "subdomain".
    """

    def __init__(self, symbol: str, parent_domain: 'DomainTerm', excluded_constants=None) -> None:
        """Excluded constants are optional but are where we specify which elements of the parent domain
        cannot appear in this domain"""
        self.symbol = symbol
        # DEBUG: making sure we don't pass in a SetOfConstants, which can't be a parent
        assert(isinstance(parent_domain, (RootDomain, DomainVariable)))
        self.parent_domain = parent_domain
        self.excluded_constants = excluded_constants if excluded_constants is not None else set()

    def difference(self, other: 'DomainTerm') -> 'DomainTerm':
        """Get the difference between this domain variable and another domain TERM"""
        raise NotImplementedError()

    @property
    def size(self) -> int:
        """Return the maximum size of this DomainVariable
        NOTE TODO: At the moment I take this as the size of the parent minus the size of the excluded
        constants (which are assumed to be allowed in the parent), but maybe something else is better"""
        return self.parent_domain.size - len(self.excluded_constants)

    def __eq__(self, other: Any) -> bool:
        """For now, two DomainVariables are equal when they have the same parent domain
        and excluded constants
        TODO: Work this out in more detail"""
        if not isinstance(other, DomainVariable):
            return False
        return self.parent_domain == other.parent_domain and self.excluded_constants == other.excluded_constants

    def __hash__(self) -> int:
        return hash((self.parent_domain, self.excluded_constants))

    def __str__(self) -> str:
        return f'{self.symbol}'

    def __repr__(self) -> str:
        return self.__str__()

