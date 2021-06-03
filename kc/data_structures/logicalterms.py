"""
Classes for types of logical terms, which are constants and variables
"""

from abc import ABC, abstractmethod

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



class LogicalVariable(LogicalTerm):
    """
    A FOL logical variable.
    """
    def __init__(self, symbol: str) -> None:
        self.symbol = symbol

    def __str__(self) -> str:
        return f'{self.symbol}'

    def __repr__(self) -> str:
        return self.__str__()


if __name__ == '__main__':
    c = Constant('a')
    print(c.value)
    print(c)

    x = LogicalVariable('X')
    print(x.symbol)
    print(x)
    print([x])



