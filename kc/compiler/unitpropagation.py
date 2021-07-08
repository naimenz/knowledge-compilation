"""File for unit propagation compilation rule"""

from kc.data_structures import AndNode, ConstrainedClause, UnconstrainedClause, ConstraintSet, UnitClause, Literal, SetOfConstants, CNF
from kc.compiler import KCRule

from typing import Optional, Tuple, List, Any, Set
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler
    from kc.data_structures import CNF, Clause, ConstrainedAtom, LogicalVariable

class UnitPropagation(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['Clause']]:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)
        Returns True and the unit clause if applicable, and False, None otherwise"""
        if len(cnf.clauses) == 1:
            return False, None # if there is only one clause, unit propagation does no good
        unit_clauses = [clause for clause in cnf.clauses if len(clause.literals) == 1]
        if len(unit_clauses) > 0:
            # we have to make sure the unit_clause actually has type UnitClause
            return True, unit_clauses[0].to_unit_clause()
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', unit_clause: 'UnitClause', compiler: 'Compiler') -> 'AndNode':
        """Apply UnitPropagation (which requires applying splitting and conditioning)
        and return an NNFNode"""
        unitpropagated_clauses: List['Clause'] = []
        u_atom = unit_clause.get_constrained_atoms()[0] # only one literal so we can access it directly
        for gamma in cnf.clauses:
            split_gammas = cls.split(gamma, u_atom)
            for gamma_s in split_gammas:
                conditioned_clause = cls.condition(gamma_s, unit_clause)
                if not conditioned_clause is None:
                    unitpropagated_clauses.append(conditioned_clause)
        propagated_cnf = CNF(unitpropagated_clauses)
        unit_cnf = CNF([unit_clause])
        return AndNode(compiler.compile(propagated_cnf), compiler.compile(unit_cnf))

    @classmethod
    def atoms_need_splitting(cls, this_atom: 'ConstrainedAtom', other_atom: 'ConstrainedAtom') -> bool:
        """Does this_atom need splitting with respect to other_atom?
         Returns True if it does, and False otherwise."""
        mgu_eq_classes = this_atom.get_constrained_atom_mgu_eq_classes(other_atom)
        independent = True if mgu_eq_classes is None else False
        return not independent and other_atom.does_not_subsume(this_atom, mgu_eq_classes)

    @classmethod
    def split(cls, gamma: 'Clause', A: 'ConstrainedAtom') -> List['Clause']:
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
                # NOTE DEBUG: before splitting, try propagating equality constraints
                gamma_mgu = ConstrainedClause(gamma.literals, joint_variables, cs_mgu).propagate_equality_constraints()
                # NOTE DEBUG: checking satisfiability first
                if gamma_mgu.cs.is_satisfiable():
                    return_clauses += cls.split(gamma_mgu, A)
            else:
                gamma_mgu = UnconstrainedClause(gamma.literals)
                return_clauses += cls.split(gamma_mgu, A)

        # loop over all constraints to negate for gamma_rest
        for e in cs_theta.join(cs_A):
            not_e = ~e
            cs_rest = cs_gamma.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            gamma_rest: 'Clause'
            if cs_rest.is_satisfiable():
                if joint_variables or cs_rest.is_non_empty():
                    # NOTE DEBUG: before splitting, try propagating equality constraints
                    gamma_rest = ConstrainedClause(gamma.literals, joint_variables, cs_rest).propagate_equality_constraints()
                    # NOTE DEBUG: checking satisfiability first
                    # NOTE: splitting the gamma_rests recursively as we build them
                    if gamma_rest.cs.is_satisfiable():
                        return_clauses += cls.split(gamma_rest, A)
                else:
                    gamma_rest = UnconstrainedClause(gamma.literals)
                    # NOTE: splitting the gamma_rests recursively as we build them
                    return_clauses += cls.split(gamma_rest, A)

        return return_clauses

    @classmethod
    def condition(cls, gamma: 'Clause', c_literal: 'UnitClause') -> Optional['Clause']:
        """Condition the constrained clause 'gamma' with respect to the unit clause (constrained literal) 'literal'.
        NOTE: assumes that gamma is split wrt the atom in 'literal'
        TODO: currently assuming there are no free or domain variables present."""
        if gamma.is_subsumed_by_literal(c_literal): # gamma is redundant when we have 'literal'
            return None
        else:
            necessary_literals = cls._discard_unsatisfied_literals(gamma, c_literal)
            Lambda: 'Clause'
            if isinstance(gamma, ConstrainedClause):
                Lambda = ConstrainedClause(necessary_literals, gamma.bound_vars, gamma.cs)
            else:
                Lambda = UnconstrainedClause(necessary_literals)
            return Lambda


    @classmethod
    def _discard_unsatisfied_literals(cls, gamma: 'Clause', c_literal: 'UnitClause') -> List['Literal']:
        """Return only the literals that are not unsatisfied by 'literal'"""
        lambdas = gamma.get_constrained_literals()
        necessary_literals = []
        for lam in lambdas:
            lam_literal = lam.literal

            if isinstance(gamma, ConstrainedClause):
                negated_constrained_literal = UnitClause([~lam_literal], gamma.bound_vars, gamma.cs)
            else:
                empty_bound_vars: Set['LogicalVariable'] = set()
                empty_cs = ConstraintSet([])
                negated_constrained_literal = UnitClause([~lam_literal], empty_bound_vars, empty_cs)

            if not negated_constrained_literal.is_subsumed_by_literal(c_literal):
                necessary_literals.append(lam_literal)
        return necessary_literals
