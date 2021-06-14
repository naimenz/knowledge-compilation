"""File for messing around with closing substitution ideas"""

from kc.data_structures import *
from kc.util import *

from typing import Set

from itertools import product

# a motivating example
X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
a, b, c = Constant('a'), Constant('b'), Constant('c')
D = SetOfConstants([a, b])
Universe = SetOfConstants([a, b, c])
p, q = Predicate('p', 1), Predicate('q', 1)
pX = Literal(Atom(p, [X]), True)
qX = Literal(Atom(q, [X]), True)

XinD = InclusionConstraint(X, D)
YeqZ = EqualityConstraint(Y, Z)
cs_gamma = ConstraintSet([XinD, ~YeqZ])
cs_A = ConstraintSet([XinD])

# sub = Substitution([(Z, Y)])
# new_cs_gamma = cs_gamma.apply_substitution(sub)
# print("hi",new_cs_gamma)

gamma = ConstrainedClause(UnconstrainedClause([pX, qX]), [X], cs_gamma)
A = ConstrainedAtom(UnconstrainedClause([pX]), [X], cs_A)

def get_free_logical_variables_in_clause(clause: 'ConstrainedClause') -> Set['LogicalVariable']:
    """Get only the free LOGICAL variables that appear in a clause.
    Free variables can appear in the constraint set OR the unconstrained clause."""

    all_variables: Set['LogicalVariable'] = set()
    for literal in clause.unconstrained_clause.literals:
        all_variables.update(term for term in literal.atom.terms if isinstance(term, LogicalVariable))
    cs_logical_variables = get_logical_variables_from_cs(clause.cs)
    all_variables.update(cs_logical_variables)
    # return only the variables that are not quantified over
    return all_variables - clause.bound_vars


def get_closing_substitutions(clause: 'ConstrainedClause', universe: 'SetOfConstants') -> Set['Substitution']:
    """Get all satisfiable closing substitutions for a clause with free logical variables given a specific universe of constants."""
    free_logical_variables = get_free_logical_variables_in_clause(clause)
    # For a closing substitution, all free variables are substituted with constants or bound variables 
    bound_vars: Set['LogicalTerm'] = cast(Set['LogicalTerm'], clause.bound_vars) # hack for type checking
    universe_constants: Set['LogicalTerm'] =  cast(Set['LogicalTerm'], universe.constants) # hack for type checking
    substitution_targets: Set['LogicalTerm'] = bound_vars.union(universe_constants)
    # construct ALL possible closing substitutions first, then weed out the ones that are unsatisfiable
    right_hand_sides = product(substitution_targets, repeat=len(free_logical_variables))
    all_substitutions = set(Substitution(zip(free_logical_variables, rhs)) for rhs in right_hand_sides)
    valid_substitutions = set(sub for sub in all_substitutions if is_satisfiable(clause.cs.apply_substitution(sub)))
    print(len(all_substitutions), len(valid_substitutions))
    return valid_substitutions

print(gamma)
print(get_free_logical_variables_in_clause(gamma))
print(get_closing_substitutions(gamma, Universe))



