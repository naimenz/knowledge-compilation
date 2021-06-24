"""
File for equivalence classes, used when computing mgus (following what was done in Forclift).
"""

from kc.data_structures import LogicalVariable, LogicalTerm, Constant, Substitution

from abc import ABC

from functools import reduce
from typing import Sequence, Set, Iterable, List, Tuple, FrozenSet, Optional
from typing import TypeVar, cast, Generic

VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']
TEC = TypeVar('TEC', bound='EquivalenceClass') # this is a type var so I can work with eq classes or var eq classes

class EquivalenceClass:
    """A class implementing an equivalence class between logical terms in FOL-DC.
    """
    def __init__(self, members: Iterable['LogicalTerm']) -> None:
        """Instantiate the private set of members that belong to this equivalence class"""
        self._members = frozenset(members)
        self._variables = frozenset(t for t in members if isinstance(t, LogicalVariable))
        self._constants = frozenset(t for t in members if isinstance(t, LogicalVariable))

    @property
    def members(self) -> FrozenSet['LogicalTerm']:
        return self._members

    def overlaps(self, other: 'EquivalenceClass') -> bool:
        """Returns True if there is an element shared by this EquivalenceClass and the other"""
        return len(self.members.intersection(other.members)) > 0

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
        return not self.is_consistent

    def join(self, other: 'EquivalenceClass') -> 'EquivalenceClass':
        """Create a larger equivalence class by adding the members of another"""
        new_members = self.members.union(other.members)
        return EquivalenceClass(new_members)

    def get_variables_only(self) -> 'VariableEquivalenceClass':
        """Take an equivalence class and remove all constants"""
        variables = [term for term in self.members if isinstance(term, LogicalVariable)]
        return VariableEquivalenceClass(variables)

    def _get_constant_or_random_variable(self) -> 'LogicalTerm':
        """If there's a constant in this equivalence class, return that.
        Otherwise, return an arbitrary variable"""
        for element in self.members:
            if isinstance(element, Constant):
                return element
        return element

    def _partition_overlapping_disjoint_classes(self,
                                                other_eq_classes: Set['TEC']
                                                ) -> Tuple[Set['TEC'], Set['TEC']]:
        """Return two lists of equivalence classes, one for those classes that overlap with this
        class, and another for the rest"""
        overlapping: Set['TEC']
        disjoint: Set['TEC']
        overlapping, disjoint = set(), set()
        for eq_class in other_eq_classes:
            if eq_class.overlaps(self):
                overlapping.add(eq_class)
            else:
                disjoint.add(eq_class)
        return overlapping, disjoint


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


class EquivalenceClasses(Generic[TEC]):
    """This is wrapper class around a set of 'EquivalenceClass' objects.
    NOTE: It uses generic typing to work with either EquivalenceClass or VariableEquivalenceClass"""
    def __init__(self, eq_classes: Iterable['TEC']) -> None:
        self._eq_classes = set(eq_classes)

    @property
    def classes(self) -> Set[TEC]:
        return self._eq_classes

    def make_self_consistent(self) -> 'EquivalenceClasses[TEC]':
        """Take an initial set of equivalence classes and iterate them until they are
        self-consistent, i.e. no overlapping equivalence classes.
        NOTE: Could be a little more efficient if we checked for inconsistency of constants 
        earlier, but that would make it less clean.
        TODO: Decide if this makes sense as the place for this function"""
        remaining_eq_classes = self.classes.copy()
        final_eq_classes: Set[TEC] = set()
        while len(remaining_eq_classes) > 0:
            current_eq_class = remaining_eq_classes.pop()

            current_class_changed = True
            while current_class_changed:
                # split the equivalence classes into those that overlap with the current class and those that don't
                overlapping, disjoint = current_eq_class._partition_overlapping_disjoint_classes(remaining_eq_classes)
                # merge the overlapping eq classes into one big eq class
                merger = lambda eq_class1, eq_class2: eq_class1.join(eq_class2)
                current_eq_class = reduce(merger, overlapping, current_eq_class)
                remaining_eq_classes = disjoint
                current_class_changed = (len(overlapping) > 0)

            final_eq_classes.add(current_eq_class)
        # return as an instance of TEC (either 
        return self.__class__(final_eq_classes)

    def to_substitution(self) -> 'Substitution':
        """Convert a list of equivalence classes into a Substitution of variables that conveys
        the same information."""
        var_term_pairs: List[VarTermPair] = []
        # we process each equivalence class in turn, generating var-term pairs
        for eq_class in self.classes:
            # choose what to map everything to
            mapping_target = eq_class._get_constant_or_random_variable()
            # now we get term-var pairs for all terms that are NOT the target
            for term in eq_class.members:
                if term != mapping_target:
                    if isinstance(term, Constant):
                        raise ValueError('There should have been at most 1 constant')
                    var = cast('LogicalVariable', term) # hack for type-checking
                    var_term_pairs.append((var, mapping_target))
        return Substitution(var_term_pairs)

    def to_constraint_set(self) -> 'ConstraintSet':
        """Convert a sequence of equivalence classes into a substitution.
        We go via a subtitution because we already have functions for that"""
        substitution = self.to_substitution()
        cs = substitution_to_constraint_set(substitution)
        return cs

    def __iter__(self):
        return iter(self.classes)

    def __str__(self) -> str:
        return f'ECs{str(self.classes)}'

    def __repr__(self) -> str:
        return self.__str__()

if __name__ == '__main__':
    ec1 = VariableEquivalenceClass([])
    ec2 = VariableEquivalenceClass([])
    # ecs = EquivalenceClasses['Type[EquivalenceClass]']([ec1, ec2])
    ecs = EquivalenceClasses([ec1, ec2])
    def f(ec: 'VariableEquivalenceClass') -> None:
        print("hi")
    f(ecs.classes.pop())
    print(type(ecs.classes.pop()))
