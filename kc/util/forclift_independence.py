"""This file is for translating Forclift's Independence and unify algorithms
into Python so I understand how they work and can write them more cleanly"""

from kc.data_structures import *
from kc.util import get_constrained_atoms

from functools import reduce

from typing import Optional, Tuple, List, Sequence, TypeVar
from typing import cast
# defining type alias to simplify type hinting
ECList = List['EquivalenceClass']
ECSeq = Sequence['EquivalenceClass']
VECList = List['VariableEquivalenceClass']
VECSeq = Sequence['VariableEquivalenceClass']
VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']
# this is a type var so I can work with eq classes or var eq classes
TEC = TypeVar('TEC', bound='EquivalenceClass')

def tryIndependentSubtheories(cnf: 'CNF'
                              ) -> Optional[Tuple['CNF', 'CNF']]:
    """They return an NNFNode but I don't have that set up yet, so I'll just
    return the independent subtheories"""

    def partition(depClauses: List['ConstrainedClause'],
                  indepClauses: List['ConstrainedClause']
                  ) -> Tuple[List['ConstrainedClause'], List['ConstrainedClause']]:
        """They use this function to construct the independent subtheories recursively"""
        if len(indepClauses) == 0:
            return depClauses, []
        # this is done with a pattern match in scala
        else:
            if len(depClauses) > 0:
                clause, rest = depClauses[0], depClauses[1:]
                indep = [other_clause for other_clause in indepClauses if clauses_independent(clause, other_clause)]
                dep = [other_clause for other_clause in indepClauses if not clauses_independent(clause, other_clause)]
                depAll, indepAll = partition(rest + dep, indep)
                return [clause] + depAll, indepAll
            else:
                return [], indepClauses

    clauses = list(cnf.clauses)
    dep, indep = partition([clauses[0]], clauses[1:])
    if len(indep) == 0:
        return None
    else:
        return CNF(dep), CNF(indep)

def clauses_independent(clause: 'ConstrainedClause', other_clause: 'ConstrainedClause') -> bool:
    """Are these two constrained clauses independent of each other?"""
    if is_conditional_contradiction(clause) or is_conditional_contradiction(other_clause):
        return True

    all_independent = True
    for c_atom in get_constrained_atoms(clause):
        for other_c_atom in get_constrained_atoms(other_clause):
            if constrained_atoms_unify(c_atom, other_c_atom):
                all_independent = False
    return all_independent

def constrained_atoms_unify(c_atom: 'ConstrainedAtom', other_c_atom: 'ConstrainedAtom') -> bool:
    """Returns True if there is a substitution that unifies c_atom and other_c_atom, otherwise False."""
    return not get_constrained_atom_mgu_substitution(c_atom, other_c_atom) is None

def get_constrained_atom_mgu_substitution(c_atom: 'ConstrainedAtom',
                                          other_c_atom: 'ConstrainedAtom'
                                          ) -> Optional['Substitution']:
    """Get the mgu of two constrained atoms.
    This is the same as the mgu for two unconstrained atoms, with an additional check 
    to see if the constraint sets conjoined with the mgu are satisfiable."""
    unconstrained_mgu = get_unconstrained_atom_mgu_substitution(c_atom.atom, other_c_atom.atom)
    if unconstrained_mgu is None:
        return None
    cs_mgu = substitution_to_constraint_set(unconstrained_mgu)
    combined_constraint_set = c_atom.cs.join(other_c_atom.cs).join(cs_mgu)
    if is_satisfiable(combined_constraint_set):
        return unconstrained_mgu
    else:
        return None

def is_satisfiable(cs: 'ConstraintSet') -> bool:
    """Is the given constraint set satisfiable? i.e. are there any substitutions to
    its variables that do not contain contradictions?
    NOTE TODO: There are two major flaws with this function still:
        1) It cannot handle domain variables
        2) If there are more mutually unequal variables in a domain than there are elements of that domain
        it is not satisfiable but this function cannot determine that yet"""
    # First, we construct variable equivalence classes from the equality constraints 
    eq_classes = get_var_eq_classes_from_cs(cs)
    # Next, we make sure that these are consistent with the inequality constraints
    if not consistent_with_inequality_constraints(eq_classes, cs):
        return False
    # Then we construct the possible domains for each equivalence class and
    # check that they are non-empty
    if not consistent_with_set_constraints(eq_classes, cs):
        return False
    # finally, we assume that it is satisfiable
    return True


def get_var_eq_classes_from_cs(cs: 'ConstraintSet') -> 'VECSeq':
    """Return the variable equivalence classes given by the equality constraints
    of a constraint set
    NOTE: This uses the assumption that equality constraints only contain logical variables"""
    initial_eq_classes = [] 
    for constraint in cs:
        if isinstance(constraint, EqualityConstraint):
            initial_eq_classes.append(VariableEquivalenceClass(constraint.terms))

    final_eq_classes = make_eq_classes_self_consistent(initial_eq_classes)
    return final_eq_classes

def consistent_with_inequality_constraints(eq_classes: Sequence['VariableEquivalenceClass'], cs: 'ConstraintSet') -> bool:
    """If the equality classes conflict with the inequality constraints, returns False
    NOTE: assuming the inequality constraints only contain logical variables"""
    for eq_class in eq_classes:
        for constraint in cs:
            if isinstance(constraint, InequalityConstraint):
                if constraint.left_term in eq_class and constraint.left_term in eq_class: 
                    return False
    return True

def consistent_with_set_constraints(eq_classes: Sequence['VariableEquivalenceClass'], cs: 'ConstraintSet') -> bool:
    """If any equality class has an empty domain, then it is not consistent
    NOTE: For now, only works with SetOfConstants, not DomainVariable"""
    for eq_class in eq_classes:
        included_domains = []
        excluded_domains = []
        for constraint in cs:
            if isinstance(constraint, InclusionConstraint) and constraint.logical_term in eq_class:
                included_domains.append(constraint.domain_term)
            elif isinstance(constraint, NotInclusionConstraint) and constraint.logical_term in eq_class:
                excluded_domains.append(constraint.domain_term)
        # hack for type checking for now since we only work with SetOfConstants
        included_domain = cast('SetOfConstants', DomainTerm.union(*included_domains)) 
        excluded_domain = cast('SetOfConstants', DomainTerm.intersection(*excluded_domains)) 
        allowed_constants = included_domain.difference(excluded_domain)
        if allowed_constants.size == 0:
            return False
    return True


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

def get_unconstrained_atom_mgu_eq_classes(atom: 'Atom', other_atom: 'Atom') -> Optional['ECList']:
    """Compute the mgu of two unconstrained atoms using equivalence classes,

    Returns a list of equivalence classes or None, if unsuccessful.
    """
    if atom.predicate != other_atom.predicate:
        return None

    # we build up equivalence classes, starting from equivalences between terms in the same position
    # of the two atoms
    final_eq_classes: 'ECList' = []
    term_pairs = zip(atom.terms, other_atom.terms)
    # we only include equivalence classes with more than one unique element
    initial_eq_classes = [EquivalenceClass([t1, t2]) for t1, t2 in term_pairs if t1 != t2]

    final_eq_classes = make_eq_classes_self_consistent(initial_eq_classes)
    if any(eq_class.is_inconsistent for eq_class in final_eq_classes):
        return None
    else:
        return final_eq_classes

def make_eq_classes_self_consistent(initial_eq_classes: Sequence['TEC']) -> List['TEC']:
    """Take an initial set of equivalence classes and iterate them until they are self-consistent,
    i.e. no overlapping equivalence classes.
    NOTE: Could be a little more efficient if we checked for inconsistency of constants earlier,
    but that would make it less clean."""
    remaining_eq_classes = initial_eq_classes
    final_eq_classes: List['TEC'] = []
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

        final_eq_classes.append(current_eq_class)
    return final_eq_classes

def partition_overlapping_disjoint_classes(current_eq_class: 'TEC',
                                                    other_eq_classes: Sequence['TEC']
                                                    ) -> Tuple[Sequence['TEC'], Sequence['TEC']]:
    """Return two lists of equivalence classes, one for those classes that overlap with eq_class,
     and another for the rest"""
    overlapping, disjoint = [], []
    for eq_class in other_eq_classes:
        if eq_class.overlaps(current_eq_class):
            overlapping.append(eq_class)
        else:
            disjoint.append(eq_class)
    return overlapping, disjoint

def eq_classes_to_substitution(eq_classes: 'ECList') -> 'Substitution':
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


def is_conditional_contradiction(clause):
    """I think this means that the clause contains no literals and so
    its grounding is empty, meaning it is independent of everything."""
    return len(clause.unconstrained_clause.literals) == 0


if __name__=='__main__':
    X = LogicalVariable('X')
    Y = LogicalVariable('Y')
    Z = LogicalVariable('Z')

    alice = Constant('alice')
    bob = Constant('bob')
    charlie = Constant('charlie')

    friends = Predicate('friends', 2)
    dislikes = Predicate('dislikes', 2)
    likes = Predicate('likes', 2)

    friendsXY = Literal(Atom(friends, [X, Y]), True)
    dislikesXY = Literal(Atom(dislikes, [X, Y]), True)
    friendsZZ = Literal(Atom(friends, [Z, Z]), True)
    likesZZ = Literal(Atom(likes, [Z, Z]), True)

    uclause1 = UnconstrainedClause([friendsXY, dislikesXY])
    uclause2 = UnconstrainedClause([~friendsZZ, likesZZ])


    People = SetOfConstants([alice, bob, charlie])
    XinPeople = InclusionConstraint(X, People)
    YinPeople = InclusionConstraint(Y, People)
    ZinPeople = InclusionConstraint(Z, People)
    XeqY = EqualityConstraint(X, Y)

    cs1 = ConstraintSet([XinPeople, YinPeople])#, ~XeqY])
    cs2 = ConstraintSet([ZinPeople])

    clause1 = ConstrainedClause(uclause1, [X, Y], cs1)
    clause2 = ConstrainedClause(uclause2, [Z], cs2)

    # print(clause1, clause2)
    cnf = CNF([clause1, clause2])
    # print(cnf)
    print(len(tryIndependentSubtheories(cnf)))
