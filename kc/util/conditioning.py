"""File for Conditioning, the second auxiliary operation for UnitProp.

NOTE: I am currently not considering the case with free or domain variables."""

from kc.data_structures import *
from kc.util import *

from typing import Optional, List

def condition(gamma: 'ConstrainedClause', literal: 'UnitClause', universe: 'SetOfConstants') -> Optional['ConstrainedClause']:
    """Condition the constrained clause 'gamma' with respect to the unit clause (constrained literal) 'literal'.
    NOTE: assumes that gamma is split wrt the atom in 'literal'
    TODO: currently assuming there are no free or domain variables present."""
    if constrained_clauses_subsumed(literal, gamma, universe): # gamma is redundant when we have 'literal'
        return None
    else:
        necessary_literals = discard_unsatisfied_literals(gamma, literal, universe)
        Lambda = ConstrainedClause(UnconstrainedClause(necessary_literals), gamma.bound_vars, gamma.cs)
        return Lambda


def discard_unsatisfied_literals(gamma: 'ConstrainedClause', literal: 'UnitClause', universe: 'SetOfConstants') -> List['Literal']:
    """Return only the literals that are not unsatisfied by 'literal'"""
    lambdas = get_constrained_literals(gamma)
    necessary_literals = []
    for lam in lambdas:
        lam_literal = lam.literal
        negated_constrained_literal = UnitClause(UnconstrainedClause([~lam_literal]), gamma.bound_vars, gamma.cs)
        if not constrained_clauses_subsumed(literal, negated_constrained_literal, universe):
            necessary_literals.append(lam_literal)
    return necessary_literals

    


