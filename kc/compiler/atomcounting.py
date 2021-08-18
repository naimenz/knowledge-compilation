"""File for atom counting compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from collections import defaultdict

from typing import Tuple, Optional, List, Set, Dict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class AtomCounting(KCRule):

    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['ConstrainedAtom']]:
        """AtomCounting is applicable if the theory is shattered
        (which should have been checked before this rule) and there is an atom in cnf with exactly one bound
        logical variable
        Returns (True, atom) if applicable, and (False, None) otherwise."""

        candidate_c_atoms = []
        for clause in sorted(cnf.clauses):
            for c_atom in sorted(clause.get_constrained_atoms()):
                common_bound_variables = set(c_atom.atom.terms).intersection(c_atom.bound_vars)
                if len(common_bound_variables) == 1:
                    candidate_c_atoms.append(c_atom)
        if len(candidate_c_atoms) > 0:
            # c_atoms are chosen in order of number of occurrences of the predicate minus size of the domain see below
            best_candidate = cls._get_best_candidate(candidate_c_atoms)
            return True, best_candidate 
        else:
            return False, None

    @classmethod
    def _get_best_candidate(cls, candidate_c_atoms: List['ConstrainedAtom']) -> 'ConstrainedAtom':
        """Get the best candidate c_atom from a list of candidates using a heuristic.
        We use the Forclift heuristic: number of occurrences of the predicate minus size of the domain, e.g
            (\forall X, X \in {alice,bob}: fun(X) OR smokes(X))
        AND (\forall X, X \in {alice,bob,charlie}: likes_haggis(X) OR smokes(X))
        smokes has scores 2 - 3 = -1 and 2 - 2 = 0
        fun has score 1 - 2 = -1
        likes_haggis 1 - 3 = - 2
        thus smokes has the highest score with 0 and is chosen
        """
        predicate_counts: Dict['Predicate', int] = defaultdict(int)
        for c_atom in candidate_c_atoms:
            predicate_counts[c_atom.atom.predicate] += 1

        ordered_atom_score_tuples = []
        for c_atom in candidate_c_atoms:
            domain = c_atom.cs.get_domain_for_variable(get_element_of_set(c_atom.bound_vars))
            ordered_atom_score_tuples.append((c_atom, predicate_counts[c_atom.atom.predicate] - domain.size()))
        return max(ordered_atom_score_tuples, key=lambda atom_score_tuple: atom_score_tuple[1])[0]

    @classmethod
    def apply(cls, cnf: 'CNF', c_atom: 'ConstrainedAtom', compiler: 'Compiler') -> 'NNFNode':
        """Apply AtomCounting and return an NNFNode"""

        atom = c_atom.atom
        # c_atom might have the form (\forall X,Y \in D: smokes(X)) due to stemming from a larger clause
        # We just want to retrieve one variable, i.e. X in this case.
        bound_var = get_element_of_set(c_atom.bound_vars.intersection(atom.variables))

        # drop all irrelevant constraints (this is like getting cs_a from the pseudocode in Alg 16 in MSc, 15 in PhD)
        variable_cs = c_atom.cs.project(c_atom)

        domain_cs, domain_variable = cls._construct_domain_cs_from_variable_cs(bound_var, variable_cs)

        # This creates constraints for the variables to be in the domains D_\top, D_\bot, respectively.
        bound_var_in_domain_variable_constraint = InclusionConstraint(bound_var, domain_variable)
        bound_var_in_complement_constraint = InclusionConstraint(bound_var, domain_variable.complement)

        true_branch_cs = variable_cs.join(ConstraintSet([bound_var_in_domain_variable_constraint]))
        false_branch_cs = variable_cs.join(ConstraintSet([bound_var_in_complement_constraint]))

        true_branch = UnitClause([Literal(atom, polarity=True)], [bound_var], true_branch_cs)
        false_branch = UnitClause([Literal(atom, polarity=False)], [bound_var], false_branch_cs)

        new_cnf = cnf.join(CNF([true_branch, false_branch]))
        new_cnf.shattered = True  # it must still be shattered at this stage due to preconditions
        new_cnf.subdivided = True  # it must still be subdivided since we don't introduce new predicates

        return ExistsNode(compiler.compile(new_cnf), [domain_variable], domain_cs)
        
    @classmethod
    def _construct_domain_cs_from_variable_cs(cls,
                                              bound_var: 'LogicalVariable',
                                              variable_cs: 'ConstraintSet',
                                              ) -> Tuple['ConstraintSet', 'DomainVariable']:
        """Build a constraint set for a new domain variable"""
        parent_domain = variable_cs.get_domain_for_variable(bound_var)
        unequal_constants = variable_cs.unequal_constants_for_variable(bound_var)  # X \neq Alice
        # we only exclude constants that could possibly be part of the domain
        excluded_constants = unequal_constants.intersection(parent_domain.possible_constants)
        subdomain_variable = DomainVariable(parent_domain.symbol + '-top',
                                            parent_domain,
                                            excluded_constants=excluded_constants)

        # the new domain variable must be a subset of the allowed domain for the bound variable
        # e.g. if we have X \neq Alice, and we are constructing X \in D then D cannot contain Alice
        var_constant_constraints = cls._build_constraints_between_domain_var_and_constants(subdomain_variable,
                                                                                           excluded_constants)
        domain_cs = (ConstraintSet([*var_constant_constraints, SubsetConstraint(subdomain_variable, parent_domain)]))
        return domain_cs, subdomain_variable

    @classmethod
    def _build_constraints_between_domain_var_and_constants(cls,
                                                            variable: 'DomainVariable',
                                                            constants: Set['Constant']
                                                            ) -> List['NotSubsetConstraint']:
        """Return a list of constraints specifying which constants are NOT contained in a DomainVariable
        e.g. Alice \notin D_\bot.
        Note that constraints involving constants work with a SetOfConstants, i.e. X \notin set(Alice)"""
        constraints: List['NotSubsetConstraint'] = []
        for constant in constants:
            constraints.append(NotSubsetConstraint(SetOfConstants([constant]), variable))
        return constraints
