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
        if not isinstance(other, Constant):
            return False
        return self.value == other.value

    def __hash__(self) -> int:
        return hash(self.__repr__())


class LogicalVariable(LogicalTerm):
    """
    A FOL logical variable.
    """
    def __init__(self, symbol: str) -> None:
        self.symbol = symbol

    def __eq__(self, other: Any) -> bool:
        """Two logical variables are equal if they have the same symbol.

        NOTE: I'm not sure if this is the right call, since symbols are sometimes reused.
        TODO: figure out what to do about symbol reuse and remove warning."""
        if not isinstance(other, LogicalVariable):
            return False
        # this is the case I'm not sure about
        if (self.symbol == other.symbol) and not (self is other):
            print(f"WARNING: Treating LogicalVariables {self.symbol} as equal but are not the same object")
            return True
        elif (self is other):
            return True
        else:
            return False

    def __str__(self) -> str:
        return f'{self.symbol}'

    def __repr__(self) -> str:
        return self.__str__()

    def __hash__(self) -> int:
        """Hashing LogicalVariables so I can use them as dict keys.

        NOTE: this may cause collisions for variables with the same symbol, but that's kind of what we want?"""
        return hash(self.__repr__())


if __name__ == '__main__':
    c = Constant('a')
    print(c.value)
    print(c)

    x = LogicalVariable('X')
    print(x.symbol)
    print(x)
    print([x])



