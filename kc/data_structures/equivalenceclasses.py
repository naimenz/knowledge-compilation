"""
File for equivalence classes, used when computing mgus (following what was done in Forclift).
"""

from kc.data_structures.logicalterms import *

from typing import Sequence, Set, Iterable

class EquivalenceClass:
    """A class implementing an equivalence class between logical terms in FOL-DC.
    """
    def __init__(self, members: Iterable['LogicalTerm']) -> None:
        """Instantiate the private set of members that belong to this equivalence class"""
        self._members = set(members)

    @property
    def members(self) -> Set['LogicalTerm']:
        return self._members

    def overlaps(self, other: 'EquivalenceClass') -> bool:
        """Returns True if there is an element shared by this EquivalenceClass and the other"""
        return len(self.members.intersection(other.members)) > 0

    @property
    def is_consistent(self) -> bool:
        """If this equivalence class contains an obvious contradiction (i.e. two different constants)
        then it is inconsistent and we return False"""
        for term1 in self.members:
            if isinstance(term1, Constant):
                for term2 in self.members:
                    if isinstance(term2, Constant) and term1 != term2:
                        return False
        return True

    @property
    def is_inconsistent(self) -> bool:
        """Negation of is_consistent (for convenience)"""
        return not self.is_consistent

    def join(self, other: 'EquivalenceClass') -> 'EquivalenceClass':
        """Create a larger equivalence class by adding the members of another"""
        new_members = self.members.union(other.members)
        return EquivalenceClass(new_members)

    def __iter__(self):
        return iter(self.members)

    def __str__(self) -> str:
        return str(self._members)

    def __repr__(self) -> str:
        return self.__str__()


class VariableEquivalenceClass(EquivalenceClass):
    """A subclass of EquivalenceClass that can contain only logical variables"""
    def __init__(self, members: Iterable['LogicalVariable']) -> None:
        super().__init__(members)

