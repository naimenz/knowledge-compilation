"""
A file to run the splitting auxiliary compilation algorithm
"""
from kc.data_structures import *
from kc.util import *

def split(gamma: 'ConstrainedClause', a: 'ConstrainedAtom') -> List['ConstrainedClause']:
    """Split the constrained clause gamma with respect to the constrained atom a.
    Returns a sequence of constrained clauses that are split with respect to a"""
    constrained_atoms = get_constrained_atoms(gamma)
    viable_atoms = [a_gamma for a_gamma in constrained_atoms if constrained_atoms_not_independent_and_not_subsumed(a, a_gamma)]
    if len(viable_atoms) == 0: # we are done if all are independent or subsumed
        return [gamma]
    a_gamma = viable_atoms[0]

    cs_gamma = a_gamma.cs
    cs_a = a.cs
    theta = get_constrained_atom_mgu_substitution(a, a_gamma)
    if theta is None:
        raise ValueError(f"a_gamma = {a_gamma} and a = {a} are independent but shouldn't be")
    cs_theta = substitution_to_constraint_set(theta)
    cs_mgu = cs_gamma.join(cs_theta).join(cs_a)
    joint_variables = a.bound_vars.union(a_gamma.bound_vars)
    gamma_mgu = ConstrainedClause(gamma.unconstrained_clause, joint_variables, cs_mgu)

    Gamma_rests: List['ConstrainedClause'] = []
    # loop over all constraints to negate for gamma_rest
    for e in cs_theta.join(cs_a):
        not_e = ~e
        cs_rest = cs_gamma.join(ConstraintSet([not_e]))
        # before going further, check if the constraint set for this clause is even satisfiable
        if is_satisfiable(cs_rest):
            gamma_rest = ConstrainedClause(gamma.unconstrained_clause, joint_variables, cs_rest)
            # NOTE: splitting the gamma_rests recursively as we build them
            Gamma_rests += split(gamma_rest, a)
    # if gamma_mgu is not satisfiable, we just return the rest
    if is_satisfiable(gamma_mgu.cs):
        return split(gamma_mgu, a) + Gamma_rests
    else:
        return Gamma_rests



