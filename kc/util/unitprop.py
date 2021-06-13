"""File to perform unit propagation.
NOTE: as with splitting and conditioning, this is currently designed to work only without free or 
domain variables."""
from kc.data_structures import *
from kc.util import *

from typing import List

def unitprop(delta: 'CNF', u: 'UnitClause') -> 'CNF':
    """Perform unit propagation on a CNF theory delta. 
    Assumes that there is a unit clause in delta.
    TODO: make this work with free and domain vars.
    TODO: Add Compile call to return statement"""
    unitpropagated_clauses: List['ConstrainedClause'] = [u] # u is not unitpropagated but we need it anyway
    u_atom = get_constrained_atoms(u)[0] # only one literal
    for gamma in delta.clauses:
        split_gammas = split(gamma, u_atom)
        for gamma_s in split_gammas:
            conditioned_clause = condition(gamma_s, u)
            if not conditioned_clause is None:
                unitpropagated_clauses.append(conditioned_clause)
    return CNF(unitpropagated_clauses)
