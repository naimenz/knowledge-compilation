"""File for unit propagation compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import *

from typing import Optional, Tuple, List, Any
class UnitPropagation(KCRule):
    
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional['ConstrainedClause']]:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)
        Returns True and the unit clause if applicable, and False, None otherwise"""
        unit_clauses = [clause for clause in delta.clauses if len(clause.unconstrained_clause.literals) == 1]
        if len(unit_clauses) > 0:
            return True, unit_clauses[0]
        return False, None

    @classmethod
    # NOTE: We annotate compiler as 'Any' to avoid circularly importing the Compiler
    def apply(cls, delta: 'CNF', unit_clause: 'ConstrainedClause', compiler: Any) -> 'NNFNode':
        """Apply UnitPropagation (which requires applying splitting and conditioning)
        and return an NNFNode"""
        raise NotImplementedError('UnitPropagation.apply not implemented')

        unitpropagated_clauses: List['ConstrainedClause'] = []
        u_atom = get_constrained_atoms(unit_clause)[0] # only one literal so we can access it directly
        for gamma in delta.clauses:
            split_gammas = cls._split(gamma, u_atom)
            for gamma_s in split_gammas:
                conditioned_clause = cls._condition(gamma_s, unit_clause)
                if not conditioned_clause is None:
                    unitpropagated_clauses.append(conditioned_clause)
        # TODO: make circuit nodes
        propagated_cnf = CNF(unitpropagated_clauses)
        unit_cnf = CNF([unit_clause])
        return AndNode(compiler.compile(propagated_cnf), compiler.compile(unit_cnf))

    @classmethod
    def _split(cls, gamma: 'ConstrainedClause', A: 'ConstrainedAtom', universe: 'SetOfConstants') -> List['ConstrainedClause']:
        """Split the constrained clause gamma with respect to the constrained atom A.
        Returns a sequence of constrained clauses that are split with respect to a"""
        constrained_atoms = get_constrained_atoms(gamma)
        viable_atoms = [a_gamma for a_gamma in constrained_atoms if constrained_atoms_not_independent_and_not_subsumed(A, a_gamma, universe)]
        if len(viable_atoms) == 0: # we are done if all are independent or subsumed
            return [gamma]
        a_gamma = viable_atoms[0]

        cs_gamma = a_gamma.cs
        cs_A = A.cs
        theta = get_constrained_atom_mgu_substitution(A, a_gamma)
        if theta is None:
            raise ValueError(f"a_gamma = {a_gamma} and A = {A} are independent but shouldn't be")
        cs_theta = substitution_to_constraint_set(theta)
        cs_mgu = cs_gamma.join(cs_theta).join(cs_A)
        joint_variables = A.bound_vars.union(a_gamma.bound_vars)
        gamma_mgu = ConstrainedClause(gamma.unconstrained_clause, joint_variables, cs_mgu)

        Gamma_rests: List['ConstrainedClause'] = []
        # loop over all constraints to negate for gamma_rest
        for e in cs_theta.join(cs_A):
            not_e = ~e
            cs_rest = cs_gamma.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            if is_satisfiable(cs_rest):
                gamma_rest = ConstrainedClause(gamma.unconstrained_clause, joint_variables, cs_rest)
                # NOTE: splitting the gamma_rests recursively as we build them
                Gamma_rests += split(gamma_rest, A, universe)
        # if gamma_mgu is not satisfiable, we just return the rest
        if is_satisfiable(gamma_mgu.cs):
            return split(gamma_mgu, A, universe) + Gamma_rests
        else:
            return Gamma_rests





