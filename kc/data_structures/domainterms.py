
"""
Classes for types of domain terms, which are sets of constants and domain variables
"""

from abc import ABC, abstractmethod
from kc.data_structures import Constant

from itertools import zip_longest

from typing import List, Set, Any, FrozenSet, Iterable
from typing import cast, Optional

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
    def union_constants(*args: 'DomainTerm') -> FrozenSet['Constant']:
        """Get the union of all the constants from the domain terms
        For ProperDomains, this uses the 'possible_constants'"""
        if len(args) == 0:
            return frozenset()
        constant_sets = []
        for arg in args:
            if isinstance(arg, ProperDomain):
                constant_sets.append(arg.possible_constants)
            elif isinstance(arg, SetOfConstants):
                constant_sets.append(arg.constants)
            all_constants: FrozenSet['Constant'] = frozenset.union(*constant_sets)
        return all_constants

    @staticmethod
    def intersect_constants(*args: 'DomainTerm') -> FrozenSet['Constant']:
        """Get the intersection of all the constants from the domain terms.
        For ProperDomains, this uses the 'possible_constants'"""
        if len(args) == 0:
            return frozenset()
        constant_sets = []
        for arg in args:
            if isinstance(arg, ProperDomain):
                constant_sets.append(arg.possible_constants)
            elif isinstance(arg, SetOfConstants):
                constant_sets.append(arg.constants)
            shared_constants: FrozenSet['Constant'] = frozenset.intersection(*constant_sets)
        return shared_constants

    # def is_subset_of(self, other: 'DomainTerm') -> bool:
    #     """Return True if THIS SetOfConstants (self) is a subset of 'other' and False otherwise
    #     NOTE: For now only works with SetOfConstants"""
    #     if all(isinstance(term, SetOfConstants) for term in (self, other) ):
    #         domains = cast(List['SetOfConstants'], (self, other) ) # hack for type checking
    #         constants: Set['Constant'] = set.intersection(*[set(domain.constants) for domain in domains])
    #         return SetOfConstants(constants) == self
    #     else:
    #         raise NotImplementedError('is_subset_of only works for SetOfConstants for now')

    # def is_superset_of(self, other: 'DomainTerm') -> bool:
    #     """Return True if THIS SetOfConstants (self) is a superset of 'other' and False otherwise
    #     NOTE: For now only works with SetOfConstants"""
    #     if all(isinstance(term, SetOfConstants) for term in (self, other) ):
    #         domains = cast(List['SetOfConstants'], (self, other) ) # hack for type checking
    #         constants: Set['Constant'] = set.intersection(*[set(domain.constants) for domain in domains])
    #         return SetOfConstants(constants) == other
    #     else:
    #         raise NotImplementedError('is_superset_of only works for SetOfConstants for now')
        

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

class ProperDomain(DomainTerm):
    """This is supposed to be an abstract class representing domain terms
    that can be the proper domain of a logical variable rather than just some constants.
    Each bound variable should have at least one ProperDomain that it belongs to.
    TODO: Work out if this approach makes sense"""

    def __init__(self, symbol: str, parent_domain: Optional['ProperDomain'], complement: 'ProperDomain'=None):
        self.symbol = symbol
        self.parent_domain = parent_domain
        self.children: List['ProperDomain'] = []  # will be updated as children are created
        self.ancestors = self.get_ancestors()
        # set as a child of the parent
        if self.parent_domain is not None:
            self.parent_domain.children.append(self)

    def get_ancestors(self) -> List['ProperDomain']:
        """Get a list of all parents of this domain, in order from most distant to most recent"""
        ancestors = []
        current_domain = self
        while current_domain.parent_domain is not None:
            current_domain = current_domain.parent_domain
            ancestors.append(current_domain)
        return list(reversed(ancestors))

    @property
    @abstractmethod
    def possible_constants(self) ->  FrozenSet['Constant']:
        """All the constants that could possibly part of this domain"""
        pass

    def intersect_with(self, other: 'ProperDomain') -> 'ProperDomain':
        """Take the intersection of two ProperDomains, which is either another ProperDomain
        or empty (None)
        NOTE TODO: anything other than empty or a ProperDomain is a 'complex intersection' which we will not
        handle because Forclift can't either"""
        # handle the simple cases first
        if isinstance(self, EmptyDomain) or isinstance(other, EmptyDomain):
            return EmptyDomain(f'{self} or {other} was empty')
        if self == other: 
            return self
        elif self.is_strict_subset_of(other):
            return self
        elif self.is_strict_superset_of(other):
            return other
        elif isinstance(self, DomainVariable) and self.complement == other:
            return EmptyDomain(f'{self} and {other} are complements')
        else:
            both_ancestors = zip_longest(self.ancestors + [self], other.ancestors + [other])
            # find the first place they differ
            for self_ancestor, other_ancestor in both_ancestors:
                if self_ancestor != other_ancestor:
                    break
            if isinstance(self_ancestor, RootDomain) or isinstance(other_ancestor, RootDomain):
                # if they differ at the root domain, then the intersection is empty (since we assumed roots are disjoint)
                return EmptyDomain(f'{self}, {other} differ at root')
        raise ValueError(f'No simple intersection found for {self}, {other}')

    @staticmethod
    def intersect_all(*args: 'ProperDomain') -> 'ProperDomain':
        """Intersect a bunch of ProperDomains together, in a pairwise fashion."""
        if len(args) == 0:
            return EmptyDomain('No arguments given: empty intersection')
        final_domain = args[0]
        for arg in args[1:]:
            final_domain = final_domain.intersect_with(arg)
        return final_domain

    def is_strict_subset_of(self, other: 'ProperDomain') -> bool:
        """Is this domain a strict subset of the other? Since RootDomains
        are disjoint, this amounts to, is this a child of the other?"""
        return other in self.ancestors

    def is_strict_superset_of(self, other: 'ProperDomain') -> bool:
        """Is this domain a strict superset of the other? Since RootDomains
        are disjoint, this amounts to, is this a parent of the other?"""
        return self in other.ancestors

class EmptyDomain(ProperDomain):
    """A domain with NO constants. 
    If this is a variable's domain, then it must be unsatisfiable."""
    def __init__(self, debug_message: str=''):
        ProperDomain.__init__(self, 'EMPTY', None)
        self.debug_message = debug_message

    def difference(self, other):
        raise NotImplementedError('EmptyDomain has no difference')

    def possible_constants(self) -> FrozenSet['Constant']:
        return frozenset()

    def size(self) -> int:
        return 0

    def __eq__(self, other: Any) -> bool:
        """All EmptyDomains are equal, because the debug message is not important"""
        return isinstance(other, EmptyDomain)

    def __hash__(self) -> int:
        """All EmptyDomains are the same"""
        return hash('EmptyDomain')

    def __str__(self) -> str:
        return f'EmptyDomain({self.debug_message})'

    def __repr__(self) -> str:
        return self.__str__()


class RootDomain(SetOfConstants, ProperDomain):
    """This is a class to represent a RootDomain, i.e. a specific set of constants
    given as part of the user-defined input. New RootDomains will not be created
    during compilation.
    NOTE: I am going to make the same assumption as Forclift: that RootDomains do not overlap"""
    def __init__(self, constants: Iterable['Constant'], symbol: str) -> None:
        SetOfConstants.__init__(self, constants)
        ProperDomain.__init__(self, symbol, None)

    @property
    def possible_constants(self) -> FrozenSet['Constant']:
        """The 'possible_constants' for a RootDomain are just all its constants"""
        return self.constants

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

class DomainVariable(ProperDomain):
    """
    A FOL domain variable.
    NOTE: This is based quite heavily on Forclift's "subdomain".
    """

    def __init__(self,
                 symbol: str,
                 parent_domain: 'ProperDomain',
                 included_constants: Iterable['Constant']=frozenset(),
                 excluded_constants: Iterable['Constant']=frozenset(),
                 complement: Optional['DomainVariable']=None
                 ) -> None:
        """Excluded constants are where we specify which elements of the parent domain
        cannot appear in this domain.
        Included constants are those that MUST be part of the domain, mostly due to constants
        that are excluded in the complement, but possibly due to inclusion constraints/
        NOTE: We now enforce that every DV has a complement and build it if it doesn't exist.  """
        ProperDomain.__init__(self, symbol, parent_domain)
        self.parent_domain = parent_domain  # hack for type checking (since it's optional in ProperDomain)
        self.excluded_constants = frozenset(excluded_constants) 
        self.included_constants = frozenset(included_constants) 
        if complement is not None:
            self.complement = complement
        else:
            # note that we swap included and excluded constants
            complement_symbol = self.parent_domain.symbol + '\\' + self.symbol
            self.complement = DomainVariable(complement_symbol, self.parent_domain, excluded_constants, included_constants, self)


    @property
    def possible_constants(self) -> FrozenSet['Constant']:
        """The 'possible_constants' for a DomainVariable are all its parent's constants
        without the ones excluded in this domain"""
        if self.parent_domain is not None:
            return self.parent_domain.possible_constants - self.excluded_constants
        else:
            raise ValueError('DomainVariable {self} without a parent!')

    def difference(self, other: 'DomainTerm') -> 'DomainTerm':
        """Get the difference between this domain variable and another domain TERM"""
        raise NotImplementedError()

    def size(self) -> int:
        """Return the maximum size of this DomainVariable
        NOTE TODO: At the moment I take this as the size of the parent minus the size of the excluded
        constants (which are assumed to be allowed in the parent), but maybe something else is better"""
        return len(self.possible_constants)

    def __eq__(self, other: Any) -> bool:
        """For now, two DomainVariables are equal when they have the same parent domain
        and excluded constants
        TODO: Work this out in more detail"""
        if not isinstance(other, DomainVariable):
            return False
        return self.parent_domain == other.parent_domain  \
                and self.excluded_constants == other.excluded_constants \
                and self.symbol == other.symbol

    def __hash__(self) -> int:
        return hash((self.symbol, self.parent_domain, self.excluded_constants))

    def __str__(self) -> str:
        return f'{self.symbol}'

    def __repr__(self) -> str:
        return self.__str__()

