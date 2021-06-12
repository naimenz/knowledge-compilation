"""
File for computing the mgu of two constrained atoms.
This relies on the EquivalenceClass class, following what was done in Forclift.
"""
from kc.data_structures import *
from kc.util.utility import *

from functools import reduce

from typing import Optional, List, Tuple
from typing import cast

# defining type alias to simplify type hinting
ECList = List['EquivalenceClass']


def get_constrained_atom_mgu_substitution(c_atom1: 'ConstrainedAtom',
                                          c_atom2: 'ConstrainedAtom'
                                          ) -> Optional['Substitution']:
    """Get the mgu of two constrained atoms.
    This is the same as the mgu for two unconstrained atoms, with an additional check 
    to see if the constraint sets conjoined with the mgu are satisfiable."""
    unconstrained_mgu = get_unconstrained_atom_mgu_substitution(c_atom1.atom, c_atom2.atom)
    if unconstrained_mgu is None:
        return None
    cs_mgu = substitution_to_constraint_set(unconstrained_mgu)
    combined_constraint_set = c_atom1.cs.join(c_atom2.cs).join(cs_mgu)
    print(combined_constraint_set)
    if is_satisfiable(combined_constraint_set):
        return unconstrained_mgu
    else:
        return None


def is_satisfiable(constraint_set: 'ConstraintSet') -> bool:
    """Check if a constraint set is satisfiable by constructing solutions to it
    and seeing if there are any"""
    # extract just the variables from each constraint
    logical_variables: Set['LogicalVariable'] = set()
    for constraint in constraint_set:
        # messy logic for isolating logical variables TODO: clean this somehow
        if isinstance(constraint, LogicalConstraint):
            if isinstance(constraint.left_term, LogicalVariable):
                logical_variables.add(constraint.left_term)
            if isinstance(constraint.right_term, LogicalVariable):
                logical_variables.add(constraint.right_term)
        elif isinstance(constraint, SetConstraint):
            if isinstance(constraint.logical_term, LogicalVariable):
                logical_variables.add(constraint.logical_term)
    # now we find solutions for those variables and see if there are any
    solutions = get_solutions(constraint_set, list(logical_variables))
    return len(solutions) > 0


def get_unconstrained_atom_mgu_substitution(atom1: 'Atom', atom2: 'Atom') -> Optional['Substitution']:
    """Compute the mgu of two unconstrained atoms using equivalence classes,
    as done in Forclift. Before returning, I convert the equivalence classes into a
    substitution.

    Returns None if the mgu doesn't exist (the atoms are independent), and a Substitution if it does.
    """
    eq_classes = get_unconstrained_atom_mgu_eq_classes(atom1, atom2)
    if not eq_classes is None:
        return eq_classes_to_substitution(eq_classes)
    else:
        return None


def get_unconstrained_atom_mgu_eq_classes(atom1: 'Atom', atom2: 'Atom') -> Optional[ECList]:
    """Compute the mgu of two unconstrained atoms using equivalence classes,

    Returns a list of equivalence classes or None, if unsuccessful.
    """
    if atom1.predicate != atom2.predicate:
        return None

    # we build up equivalence classes, starting from equivalences between terms in the same position
    # of the two atoms
    final_eq_classes: ECList = []
    term_pairs = zip(atom1.terms, atom2.terms)
    # we only include equivalence classes with more than one unique element
    initial_eq_classes = [EquivalenceClass([t1, t2]) for t1, t2 in term_pairs if t1 != t2]

    if any(eq_class.is_inconsistent for eq_class in initial_eq_classes):
        return None

    remaining_eq_classes = initial_eq_classes
    while len(remaining_eq_classes) > 0:
        current_eq_class = remaining_eq_classes[0]
        remaining_eq_classes = remaining_eq_classes[1:]

        current_class_changed = True
        while current_class_changed:
            # split the equivalence classes into those that overlap with the current class and those that don't
            overlapping, disjoint = partition_overlapping_disjoint_classes(current_eq_class, remaining_eq_classes)
            # merge the overlapping eq classes into one big eq class
            merger = lambda eq_class1, eq_class2: eq_class1.join(eq_class2)
            current_eq_class = reduce(merger, overlapping, current_eq_class)
            remaining_eq_classes = disjoint
            current_class_changed = (len(overlapping) > 0)

            if current_eq_class.is_inconsistent:
                return None

        final_eq_classes.append(current_eq_class)
    return final_eq_classes


def eq_classes_to_substitution(eq_classes: ECList) -> 'Substitution':
    """Convert a list of equivalence classes into a Substitution of variables that conveys
    the same information."""
    var_term_pairs: List[VarTermPair] = []
    # we process each equivalence class in turn, generating var-term pairs
    for eq_class in eq_classes:
        # choose what to map everything to
        mapping_target = get_constant_or_random_variable(eq_class)
        # now we get term-var pairs for all terms that are NOT the target
        for term in eq_class.members:
            if term != mapping_target:
                if isinstance(term, Constant):
                    raise ValueError('There should have been at most 1 constant')
                var = cast('LogicalVariable', term) # hack for type-checking
                var_term_pairs.append((var, mapping_target))
    return Substitution(var_term_pairs)


def substitution_to_constraint_set(substitution: 'Substitution') -> 'ConstraintSet':
    """Convert a substitution to an equivalent constraint set.

    This involves simply defining a constraint for each mapping."""
    constraints: List['Constraint'] = []
    for mapping in substitution:
        constraint = EqualityConstraint(mapping[0], mapping[1])
        constraints.append(constraint)
    return ConstraintSet(constraints)
    

def get_constant_or_random_variable(eq_class: 'EquivalenceClass') -> 'LogicalTerm':
    """If there's a constant in the equivalence class, return that.
    Otherwise, return an arbitrary variable"""
    for element in eq_class.members:
        if isinstance(element, Constant):
            return element
    return element


def partition_overlapping_disjoint_classes(current_eq_class: 'EquivalenceClass',
                                                    other_eq_classes: ECList
                                                    ) -> Tuple[ECList, ECList]:
    """Return two lists of equivalence classes, one for those classes that overlap with eq_class,
     and another for the rest"""
    overlapping, disjoint = [], []
    for eq_class in other_eq_classes:
        if eq_class.overlaps(current_eq_class):
            overlapping.append(eq_class)
        else:
            disjoint.append(eq_class)
    return overlapping, disjoint

