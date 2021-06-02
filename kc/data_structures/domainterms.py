
"""
Classes for types of domain terms, which are sets of constants and domain variables
"""

from abc import ABC, abstractmethod
from kc.data_structures.logicalterms import *

from typing import List

class DomainTerm(ABC):
    """
    An abstract base class for domain terms.
    Terms are either sets of constants or domain variables, so can never be instantiated directly.
    """


class SetOfConstants(DomainTerm):
    """
    A set of FOL constants. 
    """
    def __init__(self, constants: List['Constant']) -> None:
        self._constants = constants

    @property
    def constants(self) -> List['Constant']:
        """We use a property because the sets of constants should
        never be changed once set."""
        return self._constants

    def __str__(self) -> str:
        constant_strs = [str(constant) for constant in self.constants]
        return f"{{{', '.join(constant_strs)}}}"


class DomainVariable(DomainTerm):
    """
    A FOL domain variable.
    """
    def __init__(self, symbol: str) -> None:
        self.symbol = symbol

    def __str__(self) -> str:
        return f'{self.symbol}'


if __name__ == '__main__':
    c1 = Constant('a')
    c2 = Constant('b')
    c3 = Constant('c')

    v1 = LogicalVariable('X')

    D = SetOfConstants([c1, c2, c3])
    print(D)
