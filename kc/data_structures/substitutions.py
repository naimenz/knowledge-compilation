"""
Classes for substitutions.
"""
from kc.data_structures.logicalterms import *

from typing import List, Tuple, Dict, Any, Iterable

# defining type alias to simplify type hinting
VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']

class Substitution:
    """A FOL substitution.
    This contains a dictionary of logical variables and their substitutions (which are terms)
    """
    def __init__(self, variable_term_pairs: List[VarTermPair]) -> None:
        """The substitution dict is private and shouldn't be changed after creation"""
        self._substitution_dict = {var: term for var, term in variable_term_pairs}

    def __getitem__(self, key: 'LogicalVariable') -> 'LogicalTerm':
        """Return the term associated with a given variable"""
        return self._substitution_dict[key]

    def mappings(self) -> Iterable[VarTermPair]:
        """Return an iterator of the mappings in this substitution."""
        return self._substitution_dict.items()

    def __iter__(self):
        return iter(self.mappings())

    def __eq__(self, other: Any) -> bool:
        """Two substitutions are equal if they have all the same (variable, term) pairs"""
        if not isinstance(other, Substitution):
            return False
        return self._substitution_dict == other._substitution_dict

    def __hash__(self) -> int:
        """Hashing so I can remove duplicates with sets"""
        return hash(tuple((key, val) for key, val in self._substitution_dict.items()))

    def __str__(self) -> str:
        rightarrow_string = ' \u2192 '
        substitution_strings = [str(v) + rightarrow_string + str(t) for v, t in self._substitution_dict.items()]
        return f"{{{', '.join(substitution_strings)}}}"

    def __repr__(self) -> str:
        return self.__str__()


