"""
Classes for substitutions.
NOTE: For now, I consider substitutions of constants to variables, not free variables to variables
"""
from kc.data_structures.logicalterms import *

from typing import List, Tuple, Dict 

class Substitution:
    """A FOL substitution.
    This contains a dictionary of logical variables and their constants
    """
    def __init__(self, variable_constant_pairs: List[Tuple['LogicalVariable', 'Constant']]) -> None:
        """The substitution dict is private and shouldn't be changed after creation"""
        self._substitution_dict = {var: const for var, const in variable_constant_pairs}

    def __getitem__(self, key: 'LogicalVariable') -> 'Constant':
        """Return the constant associated with a given variable"""
        return self._substitution_dict[key]

    def __str__(self) -> str:
        rightarrow_string = ' \u2192 '
        substitution_strings = [str(v) + rightarrow_string + str(c) for v, c in self._substitution_dict.items()]
        return f"{{{', '.join(substitution_strings)}}}"

    def __repr__(self) -> str:
        return self.__str__()


if __name__ == '__main__': 
    variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]
    constants = [Constant('a'), Constant('b'), Constant('a')]

    pairs = [(v, c) for v, c in zip(variables, constants)]
    substitution = Substitution(pairs)
    print(substitution)
    print(substitution[variables[0]])
    



