"""File for unit propagation compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Optional, Tuple, List, Any
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class UnitPropagation(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['ConstrainedClause']]:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)
        Returns True and the unit clause if applicable, and False, None otherwise"""
        unit_clauses = [clause for clause in cnf.clauses if len(clause.unconstrained_clause.literals) == 1]
        if len(unit_clauses) > 0:
            return True, unit_clauses[0]
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', unit_clause: 'ConstrainedClause', compiler: 'Compiler') -> 'NNFNode':
        """Apply UnitPropagation (which requires applying splitting and conditioning)
        and return an NNFNode"""
        raise NotImplementedError('UnitPropagation.apply not implemented')

        unitpropagated_clauses: List['ConstrainedClause'] = []
        u_atom = unit_clause.get_constrained_atoms()[0] # only one literal so we can access it directly
        for gamma in cnf.clauses:
            split_gammas = cls._split(gamma, u_atom)
            for gamma_s in split_gammas:
                conditioned_clause = cls._condition(gamma_s, unit_clause)
                if not conditioned_clause is None:
                    unitpropagated_clauses.append(conditioned_clause)
        propagated_cnf = CNF(unitpropagated_clauses)
        unit_cnf = CNF([unit_clause])
        return AndNode(compiler.compile(propagated_cnf), compiler.compile(unit_cnf))

    @classmethod
    def _split(cls, gamma: 'ConstrainedClause', A: 'ConstrainedAtom', universe: 'SetOfConstants') -> List['ConstrainedClause']:
        """Split the constrained clause gamma with respect to the constrained atom A.
        Returns a sequence of constrained clauses that are split with respect to a"""
        constrained_atoms = gamma.get_constrained_atoms()
        viable_atoms = [a_gamma for a_gamma in constrained_atoms if constrained_atoms_not_independent_and_not_subsumed(A, a_gamma, universe)]
        if len(viable_atoms) == 0: # we are done if all are independent or subsumed
            return [gamma]
        a_gamma = viable_atoms[0]

        cs_gamma = a_gamma.cs
        cs_A = A.cs
        theta = A.get_constrained_atom_mgu_substitution(a_gamma)
        if theta is None:
            raise ValueError(f"a_gamma = {a_gamma} and A = {A} are independent but shouldn't be")
        cs_theta = theta.to_constraint_set()
        cs_mgu = cs_gamma.join(cs_theta).join(cs_A)
        joint_variables = A.bound_vars.union(a_gamma.bound_vars)
        gamma_mgu = ConstrainedClause(gamma.unconstrained_clause, joint_variables, cs_mgu)

        Gamma_rests: List['ConstrainedClause'] = []
        # loop over all constraints to negate for gamma_rest
        for e in cs_theta.join(cs_A):
            not_e = ~e
            cs_rest = cs_gamma.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            if cs_rest.is_satisfiable():
                gamma_rest = ConstrainedClause(gamma.unconstrained_clause, joint_variables, cs_rest)
                # NOTE: splitting the gamma_rests recursively as we build them
                Gamma_rests += cls._split(gamma_rest, A, universe)
        # if gamma_mgu is not satisfiable, we just return the rest
        if gamma_mgu.cs.is_satisfiable():
            return cls._split(gamma_mgu, A, universe) + Gamma_rests
        else:
            return Gamma_rests





