"""
Classes for substitutions.
"""
from kc.data_structures.logicalterms import *

from typing import List, Tuple, Dict, Any

class Substitution:
    """A FOL substitution.
    This contains a dictionary of logical variables and their substitutions (which are terms)
    """
    def __init__(self, variable_term_pairs: List[Tuple['LogicalVariable', 'LogicalTerm']]) -> None:
        """The substitution dict is private and shouldn't be changed after creation"""
        self._substitution_dict = {var: term for var, term in variable_term_pairs}

    def __getitem__(self, key: 'LogicalVariable') -> 'LogicalTerm':
        """Return the term associated with a given variable"""
        return self._substitution_dict[key]

    def __eq__(self, other: Any) -> bool:
        """Two substitutions are equal if they have all the same (variable, term) pairs"""
        if not isinstance(other, Substitution):
            return False
        return self._substitution_dict == other._substitution_dict
    def __str__(self) -> str:
        rightarrow_string = ' \u2192 '
        substitution_strings = [str(v) + rightarrow_string + str(t) for v, t in self._substitution_dict.items()]
        return f"{{{', '.join(substitution_strings)}}}"

    def __repr__(self) -> str:
        return self.__str__()


if __name__ == '__main__': 
    variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]
    sub_constants: List['LogicalTerm'] = [Constant('a'), Constant('b'), Constant('a')]

    pairs = [(v, c) for v, c in zip(variables, sub_constants)]
    substitution = Substitution(pairs)
    print(substitution)
    print(substitution[variables[0]])
    print(substitution == substitution)
    sub2 = Substitution(pairs[1::-1])
    print(sub2)
    print(substitution == sub2)
    



