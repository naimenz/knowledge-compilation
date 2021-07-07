"""File for atom counting compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule, ShatteredCompilation
from kc.util import get_element_of_set

from typing import Tuple, Optional
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
        for clause in cnf.clauses:
            for c_atom in clause.get_constrained_atoms():
                overlap = set(c_atom.atom.terms).intersection(c_atom.bound_vars) 
                if len(overlap) == 1:
                    return True, c_atom
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', c_atom: 'ConstrainedAtom', compiler: 'Compiler') -> 'NNFNode':
        """Apply AtomCounting and return an NNFNode"""
        atom = c_atom.atom
        bound_var = get_element_of_set(c_atom.bound_vars)  # there's only one by the preconditions
        constants = cnf.get_constants()
        free_vars = cnf.get_free_logical_variables()
        domain_terms = cnf.get_domain_terms()

        shattervar_constraint_sets = ShatteredCompilation.shatter_var(bound_var, (free_vars, constants), domain_terms)
        
        variable_cs: Optional['ConstraintSet'] = None  # this is cs_a in the pseudocode

        for cs in shattervar_constraint_sets:
            if c_atom.is_subsumed_by_c_atom( ConstrainedAtom(c_atom.literals, [bound_var], cs) ):
                variable_cs = cs

        # DEBUG
        if variable_cs is None:
            raise ValueError("No shattered cs subsumes the c_atom")

        domain_cs, domain_variable = cls._construct_domain_cs_from_variable_cs(cnf, variable_cs, bound_var)
        bound_var_in_domain_variable = InclusionConstraint(bound_var, domain_variable)

        true_branch_cs = variable_cs.join(ConstraintSet([bound_var_in_domain_variable]))
        false_branch_cs = variable_cs.join(ConstraintSet([~bound_var_in_domain_variable]))

        true_branch = UnitClause([Literal(atom, polarity=True)], [bound_var], true_branch_cs)
        false_branch = UnitClause([Literal(atom, polarity=False)], [bound_var], false_branch_cs)

        new_cnf = cnf.join(CNF([true_branch, false_branch]))
        new_cnf.shattered = True  # it must still be shattered at this stage due to preconditions

        return ExistsNode(compiler.compile(new_cnf), [domain_variable], domain_cs)
        
    @classmethod
    def _construct_domain_cs_from_variable_cs(cls,
            cnf: 'CNF',
            variable_cs: 'ConstraintSet',
            variable: 'LogicalVariable',
            ) -> Tuple['ConstraintSet', 'DomainVariable']: 
        """Build a constraint set for a new domain variable"""
        single_variable_eq_class = EquivalenceClass([variable])
        parent_domain = single_variable_eq_class.get_shared_domain_from_cs(variable_cs)
        subdomain_variable = cnf.get_new_domain_variable('D', parent_domain)

        # TODO: work out the constraints here
        domain_cs = ConstraintSet([])
        # for constraint in variable_cs:
        return domain_cs, subdomain_variable




    




