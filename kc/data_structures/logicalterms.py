"""
Classes for types of logical terms, which are constants and variables
"""

from abc import ABC, abstractmethod
from typing import Any

class LogicalTerm(ABC):
    """
    An abstract base class for logical terms.
    Terms are either constants or variables, so can never be instantiated directly.
    """
    @abstractmethod
    def __lt__(self, other: Any) -> bool:
        """All logical terms need to be comparable"""


class Constant(LogicalTerm):
    """
    A FOL constant. 
    """
    def __init__(self, value: str) -> None:
        self._value = value

    @property
    def value(self) -> str:
        """We use a property because the value of a constant should
        never be changed once set."""
        return self._value

    def __str__(self) -> str:
        return f'{self.value}'
    
    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other: Any) -> bool:
        """We say constants are equal if they have the same value"""
        return isinstance(other, Constant) and self.value == other.value

    def __hash__(self) -> int:
        return hash(self.__repr__())

    def __lt__(self, other: Any) -> bool:
        """Ordering is arbitrary - just use ordering on values"""
        if isinstance(other, Constant):
            return self.value < other.value
        elif isinstance(other, LogicalVariable):
            return self.value < other.symbol
        else:
            raise NotImplementedError(f'Cannot compare Constant and {type(other)}')


class LogicalVariable(LogicalTerm):
    """
    A FOL logical variable.
    """
    def __init__(self, symbol: str) -> None:
        self.symbol = symbol

    def __eq__(self, other: Any) -> bool:
        """Two logical variables are equal if they have the same symbol."""
        return isinstance(other, LogicalVariable) and self.symbol == other.symbol

    def __str__(self) -> str:
        return self.symbol

    def __repr__(self) -> str:
        return self.__str__()

    def __hash__(self) -> int:
        """Hashing LogicalVariables so I can use them as dict keys and in sets."""
        return hash(("LogicalVariable", self.symbol))

    def __lt__(self, other: Any) -> bool:
        """Ordering is arbitrary - just use ordering on symbols"""
        if isinstance(other, LogicalVariable):
            return self.symbol < other.symbol
        elif isinstance(other, Constant):
            return self.symbol < other.value
        else:
            raise NotImplementedError(f'Cannot compare LogicalVariable and {type(other)}')


class FreeVariable(LogicalVariable):
    """NOTE: Experimental subclass of variables specifically for FreeVariables.
    For now, doesn't have much functionality, just so we can check if something is Free
    """

    def __eq__(self, other: Any) -> bool:
        """Two free variables are equal if they have the same symbol."""
        return isinstance(other, FreeVariable) and self.symbol == other.symbol

    def __hash__(self) -> int:
        """Hashing FreeVariables so I can use them as dict keys and in sets."""
        return hash(("FreeVariable", self.symbol))
