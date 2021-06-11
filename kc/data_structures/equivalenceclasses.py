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

    def is_consistent(self) -> bool:
        """If this equivalence class contains an obvious contradiction (i.e. two different constants)
        then it is inconsistent and we return False"""
        for term1 in self.members:
            if isinstance(term1, Constant):
                for term2 in self.members:
                    if isinstance(term2, Constant) and term1 != term2:
                        return False
        return True

    def is_inconsistent(self) -> bool:
        """Negation of is_consistent (for convenience)"""
        return not self.is_consistent()

    def join(self, other: 'EquivalenceClass') -> 'EquivalenceClass':
        """Create a larger equivalence class by adding the members of another"""
        new_members = self.members.union(other.members)
        return EquivalenceClass(new_members)

    def __str__(self) -> str:
        return str(self._members)

    def __repr__(self) -> str:
        return self.__str__()

if __name__ == '__main__': 
    a, b, c = Constant('a'), Constant('b'), Constant('c')
    X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z') 
    eq_class1 = EquivalenceClass([a, b])
    eq_class2 = EquivalenceClass([X, Y])
    print(eq_class1.is_inconsistent())
    print(eq_class2.is_inconsistent())
    print(eq_class1.join(eq_class2).is_consistent())


