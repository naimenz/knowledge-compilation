"""File for unit propagation compilation rule"""

from kc.data_structures import AndNode, ConstrainedClause, UnconstrainedClause, ConstraintSet, UnitClause, Literal, SetOfConstants, CNF, Substitution, EquivalenceClasses, LogicalVariable
from kc.compiler import KCRule
DEBUG_FLAG = False

from typing import Optional, Tuple, List, Any, Set, FrozenSet
from typing import TypeVar
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from kc.compiler import Compiler
    from kc.data_structures import CNF, Clause, ConstrainedAtom, LogicalVariable

CC = TypeVar('CC', bound='ConstrainedClause')
CC2 = TypeVar('CC2', bound='ConstrainedClause')  # I need a second typevar for some functions

class UnitPropagation(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['Clause']]:
        """UnitPropagation is applicable if the theory contains a unit clause 
        (a clause with a single literal)
        Returns True and the unit clause if applicable, and False, None otherwise"""
        if len(cnf.clauses) == 1:
            return False, None # if there is only one clause, unit propagation does no good
        unit_clauses = sorted([clause for clause in cnf.clauses if len(clause.literals) == 1])
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
        for gamma in sorted(cnf.clauses):
            split_gammas = cls.split(gamma, u_atom)
            for gamma_s in sorted(split_gammas):
                conditioned_clause = cls.condition(gamma_s, unit_clause)
                if not conditioned_clause is None:
                    unitpropagated_clauses.append(conditioned_clause)
        propagated_cnf = CNF(unitpropagated_clauses, shattered=cnf.shattered)
        # NOTE: a quick check here to see if the unit_clause is really unconstrained, and if so, return
        # it that way
        if len(unit_clause.bound_vars) == 0 and unit_clause.cs == ConstraintSet([]):
            unit_cnf = CNF([UnconstrainedClause(unit_clause.literals)], shattered=cnf.shattered)
        else:
            unit_cnf = CNF([unit_clause], shattered=cnf.shattered)
        return AndNode(compiler.compile(propagated_cnf), compiler.compile(unit_cnf))

    @classmethod
    def split(cls, gamma: 'Clause', A: 'ConstrainedAtom') -> List['Clause']:
        """Split the constrained clause gamma with respect to the constrained atom A.
        Returns a sequence of constrained clauses that are split with respect to a"""
        constrained_atoms = gamma.get_constrained_atoms()
        viable_atoms = sorted([a_gamma for a_gamma in constrained_atoms if a_gamma.needs_splitting(A)])
        if len(viable_atoms) == 0: # we are done if all are independent or subsumed
            return [gamma]
        a_gamma = viable_atoms[0]

        # if we have a viable atom, apply some preprocessing to the clauses to 
        # avoid variable name issues
        # DEBUG TODO: This is experimental
        if isinstance(gamma, ConstrainedClause):
            a_gamma, gamma = cls._align_variables(A, a_gamma, gamma)
        cs_gamma = a_gamma.cs
        cs_A = A.cs

        mgu_eq_classes = A.get_constrained_atom_mgu_eq_classes(a_gamma)
        if mgu_eq_classes is None:
            raise ValueError(f"a_gamma = {a_gamma} and A = {A} are independent but shouldn't be")

        theta = mgu_eq_classes.to_substitution()
        cs_theta = theta.to_constraint_set()
        cs_mgu = cs_gamma.join(cs_theta).join(cs_A)
        joint_variables = A.bound_vars.union(a_gamma.bound_vars)

        return_clauses: List['Clause'] = []
        gamma_mgu: 'Clause'
        # print(f'{cs_mgu = }, {cs_mgu.is_satisfiable() = }')
        if cs_mgu.is_satisfiable():
            if joint_variables or cs_mgu.is_non_empty(): 
                gamma_mgu = ConstrainedClause(gamma.literals, joint_variables, cs_mgu).propagate_equality_constraints()
                if gamma_mgu.cs.is_satisfiable():
                    return_clauses += cls.split(gamma_mgu, A)
            else:
                gamma_mgu = UnconstrainedClause(gamma.literals)
                return_clauses += cls.split(gamma_mgu, A)

        # loop over all constraints to negate for gamma_rest
        for e in sorted(cs_theta.join(cs_A)):
            # NOTE DEBUG: Trying - only include constraint if it is relevant to the clause
            if not (e.terms[0] in gamma.all_variables or e.terms[1] in gamma.all_variables):
                continue
            not_e = ~e
            cs_rest = cs_gamma.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            gamma_rest: 'Clause'
            if cs_rest.is_satisfiable():
                if joint_variables or cs_rest.is_non_empty():
                    gamma_rest = ConstrainedClause(gamma.literals, joint_variables, cs_rest).propagate_equality_constraints()
                    if gamma_rest.cs.is_satisfiable():
                        return_clauses += cls.split(gamma_rest, A)
                else:
                    gamma_rest = UnconstrainedClause(gamma.literals)
                    return_clauses += cls.split(gamma_rest, A)
        return return_clauses

    @classmethod
    def _make_variables_different(cls,
                                clause: CC,
                                other_clause: CC,
                                target_clause: Optional[CC2]=None
                                ) -> Tuple[CC, Optional[CC2]]:
        """Apply preprocessing to the BOUND variables of two clauses to make them 
        all distinct. Also optionally apply this same substitution to another clause"""

        # all the variables that need to be substituted
        overlapping_variables: FrozenSet['LogicalVariable'] = clause.bound_vars.intersection(other_clause.bound_vars)

        old_other_clause = other_clause
        for variable in sorted(overlapping_variables):
            temp_cnf = CNF([clause, old_other_clause], names=None)  # taking advantage of existing methods in CNF
            sub_target = temp_cnf.get_new_logical_variable(variable.symbol)
            sub = Substitution([(variable, sub_target)])
            new_other_clause = old_other_clause.substitute(sub)
            if new_other_clause is None:
                raise ValueError('Substitution made unsatisfiable')
            else:
                old_other_clause = new_other_clause
            if target_clause is not None:
                target_clause = target_clause.substitute(sub)
        return old_other_clause, target_clause

    @classmethod
    def _align_variables(cls,
                         c_atom: 'ConstrainedAtom',
                         other_c_atom: 'ConstrainedAtom',
                         clause: 'ConstrainedClause'
                         ) -> Tuple['ConstrainedAtom', 'ConstrainedClause']:
        """'Line up' the variables in the  other_c_atom to match those in the first (c_atom), and
        apply this to the whole clause.
        This is done to make the mgu meaningful, rather than just having it rename variables to match.
        First, separate out all the bound variables.
        Then for each bound variable in the terms of this c_atom, we rename the other to match, as long as it 
        also is a bound variable."""
        other_c_atom, returned_clause = cls._make_variables_different(c_atom, other_c_atom, clause)
        assert(returned_clause is not None)  # hack for type checking
        old_clause = returned_clause
        old_other_c_atom = other_c_atom
        for term, other_term in zip(c_atom.atom.terms, other_c_atom.atom.terms):
            if term in c_atom.bound_vars and other_term in other_c_atom.bound_vars:
                assert(isinstance(other_term, LogicalVariable))  # hack for type checking
                sub = Substitution([(other_term, term)])
                new_clause = old_clause.substitute(sub)
                new_other_c_atom = old_other_c_atom.substitute(sub)
                if new_other_c_atom is None or new_clause is None:
                    raise ValueError("Shouldn't be unsatisfiable here")
                old_clause = new_clause
                old_other_c_atom = new_other_c_atom
        return old_other_c_atom, old_clause

    @classmethod
    def condition(cls, gamma: 'Clause', c_literal: 'UnitClause') -> Optional['Clause']:
        """Condition the constrained clause 'gamma' with respect to the unit clause (constrained literal) 'c_literal'.
        NOTE: assumes that gamma is split wrt the atom in 'c_literal'"""
        # print("\n==============================")
        # print(f"{gamma = }")
        # print(f"{c_literal = }")
        if gamma.is_subsumed_by_literal(c_literal): # gamma is redundant when we have 'literal'
            # print(f"SUBSUMED!")
            return None
        else:
            # print(f"NOT SUBSUMED!")
            necessary_literals = cls._discard_unsatisfied_literals(gamma, c_literal)
            Lambda: 'Clause'
            if isinstance(gamma, ConstrainedClause):
                Lambda = ConstrainedClause(necessary_literals, gamma.bound_vars, gamma.cs)
            else:
                Lambda = UnconstrainedClause(necessary_literals)
            # print(f"REDUCED TO {Lambda}!")
            return Lambda
        # print("==============================\n")


    @classmethod
    def _discard_unsatisfied_literals(cls, gamma: 'Clause', c_literal: 'UnitClause') -> List['Literal']:
        """Return only the literals that are not unsatisfied by 'literal'"""
        lambdas = gamma.get_constrained_literals()
        necessary_literals = []
        for lam in sorted(lambdas):
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
