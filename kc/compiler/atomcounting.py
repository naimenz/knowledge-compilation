"""File for atom counting compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule, ShatteredCompilation
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
        Returns True and the atom if applicable, and False, None otherwise."""

        # TODO: Heuristic for deciding which c_atom to use (from Forclift)
        # Currently using the heuristic of: choose the c_atom with the fewest terms
        candidate_c_atoms = []
        for clause in sorted(cnf.clauses):
            for c_atom in sorted(clause.get_constrained_atoms()):
                overlap = set(c_atom.atom.terms).intersection(c_atom.bound_vars) 
                if len(overlap) == 1:
                    candidate_c_atoms.append(c_atom)
        if len(candidate_c_atoms) > 0:
            best_candidate = cls._get_best_candidate(candidate_c_atoms)
            return True, best_candidate 
        else:
            return False, None

    @classmethod
    def _get_best_candidate(cls, candidate_c_atoms: List['ConstrainedAtom']) -> 'ConstrainedAtom':
        """Get the best candidate c_atom from a list of candidates using a heuristic.
        We use the Forclift heuristic: number of occurrences of the predicate minus size of the domain"""
        ordering_tuples = []
        # getting the predicate counts
        predicate_counts: Dict['Predicate', int] = defaultdict(int)
        for c_atom in candidate_c_atoms:
            predicate_counts[c_atom.atom.predicate] += 1

        for c_atom in candidate_c_atoms:
            domain = c_atom.cs.get_domain_for_variable(get_element_of_set(c_atom.bound_vars))
            ordering_tuples.append((c_atom, predicate_counts[c_atom.atom.predicate] - domain.size()))
        return max(ordering_tuples, key=lambda t: t[1])[0]

    @classmethod
    def apply(cls, cnf: 'CNF', c_atom: 'ConstrainedAtom', compiler: 'Compiler') -> 'NNFNode':
        """Apply AtomCounting and return an NNFNode
        NOTE: Trying AC in a style similar to how Forclift does"""
        atom = c_atom.atom
        bound_var = get_element_of_set(c_atom.bound_vars.intersection(atom.variables))
        constants = cnf.get_constants()
        free_vars = cnf.get_free_logical_variables()
        domain_terms = cnf.get_domain_terms()

        # drop all irrelevant constraints (this is like getting cs_a from the pseudocode in the PhD)        
        variable_cs = c_atom.cs.project(c_atom)

        domain_cs, domain_variable = cls._construct_domain_cs_from_variable_cs(cnf, bound_var, variable_cs, bound_var)
        bound_var_in_domain_variable = InclusionConstraint(bound_var, domain_variable)
        # this is equivalent to the negation of bound_var_in_domain_variable
        bound_var_in_complement = InclusionConstraint(bound_var, domain_variable.complement)

        true_branch_cs = variable_cs.join(ConstraintSet([bound_var_in_domain_variable]))
        false_branch_cs = variable_cs.join(ConstraintSet([bound_var_in_complement]))

        true_branch = UnitClause([Literal(atom, polarity=True)], [bound_var], true_branch_cs)
        false_branch = UnitClause([Literal(atom, polarity=False)], [bound_var], false_branch_cs)

        new_cnf = cnf.join(CNF([true_branch, false_branch]))
        new_cnf.shattered = True  # it must still be shattered at this stage due to preconditions

        # print("============== BEG DEBUG ===================")
        # print(f"Domain Variable after AC:\n{domain_variable}")
        # print(f"Domain CS after AC:\n{domain_cs}")
        # print(f"Theory after AC:\n{new_cnf}")
        # print("============== END DEBUG ===================")
        # raise NotImplementedError("end")
        return ExistsNode(compiler.compile(new_cnf), [domain_variable], domain_cs)
        
    @classmethod
    def _construct_domain_cs_from_variable_cs(cls,
            cnf: 'CNF',
            bound_var: 'LogicalVariable',
            variable_cs: 'ConstraintSet',
            variable: 'LogicalVariable',
            ) -> Tuple['ConstraintSet', 'DomainVariable']: 
        """Build a constraint set for a new domain variable"""
        parent_domain = variable_cs.get_domain_for_variable(bound_var)
        unequal_constants = variable_cs.unequal_constants_for_variable(bound_var)
        # we only exclude constants that could possibly be part of the domain
        excluded_constants = unequal_constants.intersection(parent_domain.possible_constants)
        # NOTE TODO: Trying out domain-top and domain-bot instead of D and underscores
        subdomain_variable = DomainVariable(parent_domain.symbol + '-top', parent_domain, excluded_constants=excluded_constants)
        # subdomain_variable = cnf.get_new_domain_variable(parent_domain.symbol, parent_domain, excluded_constants)

        # the new domain variable must be a subset of the allowed domain for the bound variable
        var_constant_constraints = cls._build_constraints_between_domain_var_and_constants(subdomain_variable, excluded_constants)
        domain_cs = (ConstraintSet([*var_constant_constraints, SubsetConstraint(subdomain_variable, parent_domain)]))
        return domain_cs, subdomain_variable

    @classmethod
    def _build_constraints_between_domain_var_and_constants(cls, variable: 'DomainVariable', constants: Set['Constant']) -> List['NotSubsetConstraint']:
        """Return a list of constraints specifying which constants are NOT contained in a DomainVariable"""
        constraints: List['NotSubsetConstraint'] = []
        for constant in constants:
            constraints.append(NotSubsetConstraint(variable, SetOfConstants([constant])))
        return constraints
