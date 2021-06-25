"""File for unit propagation compilation rule"""

from kc.data_structures import AndNode, ConstrainedClause, UnconstrainedClause, ConstraintSet
from kc.compiler import KCRule

from typing import Optional, Tuple, List, Any
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler
    from kc.data_structures import CNF, Clause, ConstrainedAtom

class UnitPropagation(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['Clause']]:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)
        Returns True and the unit clause if applicable, and False, None otherwise"""
        unit_clauses = [clause for clause in cnf.clauses if len(clause.literals) == 1]
        if len(unit_clauses) > 0:
            return True, unit_clauses[0]
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', unit_clause: 'Clause', compiler: 'Compiler') -> 'AndNode':
        """Apply UnitPropagation (which requires applying splitting and conditioning)
        and return an NNFNode"""
        unitpropagated_clauses: List['Clause'] = []
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
    def _split(cls, gamma: 'Clause', A: 'ConstrainedAtom') -> List['Clause']:
        """Split the constrained clause gamma with respect to the constrained atom A.
        Returns a sequence of constrained clauses that are split with respect to a"""
        constrained_atoms = gamma.get_constrained_atoms()
        viable_atoms = [a_gamma for a_gamma in constrained_atoms if not a_gamma.independent_or_subsumed_by(A)]
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

        return_clauses: List['Clause'] = []
        gamma_mgu: 'Clause'
        if cs_mgu.is_satisfiable():
            if joint_variables or cs_mgu.is_non_empty(): 
                gamma_mgu = ConstrainedClause(gamma.literals, joint_variables, cs_mgu)
            else:
                gamma_mgu = UnconstrainedClause(gamma.literals)
            return_clauses += cls._split(gamma_mgu, A)

        # loop over all constraints to negate for gamma_rest
        for e in cs_theta.join(cs_A):
            not_e = ~e
            cs_rest = cs_gamma.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            gamma_rest: 'Clause'
            if cs_rest.is_satisfiable():
                if joint_variables or cs_rest.is_non_empty():
                    gamma_rest = ConstrainedClause(gamma.literals, joint_variables, cs_rest)
                else:
                    gamma_rest = UnconstrainedClause(gamma.literals)
                # NOTE: splitting the gamma_rests recursively as we build them
                return_clauses += cls._split(gamma_rest, A)

        return return_clauses
