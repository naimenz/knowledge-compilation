"""File for shatter compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import powerset

from itertools import product
from functools import reduce

from typing import Tuple, Set, Iterable, Sequence
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from kc.compiler import Compiler

class ShatteredCompilation(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, None]:
        """ShatteredCompilation is applicable if the theory is not already shattered (which is represented by a flag )
        Returns True or False depending on the flag, plus None (no stored data needed)"""
        shattering_applicable = not cnf.shattered
        return shattering_applicable, None

    @classmethod
    def apply(cls, cnf: 'CNF', stored_data: None, compiler: 'Compiler') -> 'NNFNode':
        """Apply ShatteredCompilation and return an NNFNode"""
        constants = cnf.get_constants()
        free_variables = cnf.get_free_logical_variables()
        terms = (free_variables, constants)
        domains = cnf.get_domain_terms()
        shattered_clauses = [cls.shatter_clause(clause, terms, domains) for clause in cnf.clauses]

    @classmethod
    def shatter_clause(cls,
                       clause: 'Clause',
                       terms: Tuple[Set['LogicalVariable'], Set['Constant']],
                       domains: Set['DomainTerm']
                       ) -> Set['ConstrainedClause']:
        """Shatter a clause with respect to a given set of terms and domains.
        NOTE: This expects a constrained clause. I'm not really sure what to do about 
        unconstrained clauses -- hopefully any unconstrained clauses won't reach shattering."""
        if not isinstance(clause, ConstrainedClause):
            raise NotImplementedError('shatter_clause only works with ConstrainedClauses')
        # shatter each variable
        literal_variables = clause.literal_variables
        individual_shatter_var_constraints = cls._build_shatter_var_constraints(literal_variables, terms, domains)

        # set up variable inequalities
        for literal in clause.literals:
            bound_literal_variables = clause.bound_vars.intersection(literal.variables)


    @classmethod
    def _build_shatter_var_constraints(cls,
                                       variables: Iterable['LogicalVariable'],
                                       terms: Tuple[Set['LogicalVariable'], Set['Constant']],
                                       domains: Set['DomainTerm']
                                       ) -> Set['ConstraintSet']:
        """Create the shatter_var constraints (CS_A in the PhD) for a set of variables"""
        individual_shatter_var_constraints = map(lambda var: cls.shatter_var(var, terms, domains), variables)
        # make each combination of individual variable shatterings
        all_shatter_var_constraints = product(*individual_shatter_var_constraints)
        # merge each combination into a single constraint set 
        shatter_var_constraints = map(cls._merge_constraint_sets, all_shatter_var_constraints)
        return set(shatter_var_constraints)

    @classmethod
    def _merge_constraint_sets(cls, cs_iterable: Iterable['ConstraintSet']) -> 'ConstraintSet':
        """Merge a collection of constraint sets into one
        TODO: consider un-nesting this for performance"""
        return reduce(lambda cs1, cs2: cs1.join(cs2), cs_iterable)
    

    @classmethod
    def shatter_var(cls, 
                    variable: 'LogicalVariable',
                    terms: Tuple[Set['LogicalVariable'], Set['Constant']],
                    domains: Set['DomainTerm']
                    ) -> Set['ConstraintSet']:
        """Return shattered constraints for 'variable' with respect to 'terms'
        and 'domains'.
        NOTE: since we represent constant inequality with InclusionConstraints,
        many of the CS_eq terms will actually by SetConstraints."""
        equality_constraint_sets: Set['ConstraintSet'] = cls._build_equality_constraint_sets(variable, terms)
        inclusion_constraint_sets: Set['ConstraintSet'] = cls._build_inclusion_constraint_sets(variable, domains)

        product_of_cs = product(equality_constraint_sets, inclusion_constraint_sets)
        return set(cs_eq.join(cs_in) for cs_eq, cs_in in product_of_cs)

    @classmethod
    def _build_equality_constraint_sets(cls,
                                        variable: 'LogicalVariable',
                                        terms: Tuple[Set['LogicalVariable'], Set['Constant']]
                                        ) -> Set['ConstraintSet']:
        """Build the cs_eq part of the shatter_var constraint sets"""
        equality_constraint_sets: Set['ConstraintSet'] = set()

        # for building up the big 'not equal to all terms' constraint set
        not_equal_constraints: Set['Constraint'] = set()

        free_variables, constants = terms[0], terms[1]
        for free_variable in free_variables:
            var_eq_constraint = EqualityConstraint(variable, free_variable)
            equality_constraint_sets.add(ConstraintSet([var_eq_constraint]))
            not_equal_constraints.add(~var_eq_constraint)
        for constant in constants:
            const_eq_constraint = InclusionConstraint(variable, SetOfConstants([constant]))
            equality_constraint_sets.add(ConstraintSet([const_eq_constraint]))
            not_equal_constraints.add(~const_eq_constraint)
        equality_constraint_sets.add(ConstraintSet(not_equal_constraints))
        return equality_constraint_sets
        

    @classmethod
    def _build_inclusion_constraint_sets(cls,
                                        variable: 'LogicalVariable',
                                        domains: Set['DomainTerm']
                                        ) -> Set['ConstraintSet']:
        """Build the cs_in part of the shatter_var constraint sets"""
        inclusion_constraint_sets: Set['ConstraintSet'] = set()

        all_domain_subsets = powerset(domains)
        for domain_subset in all_domain_subsets:
            current_domain_constraints: Set['SetConstraint'] = set()
            for domain in domain_subset:
                variable_in_domain_constraint = InclusionConstraint(variable, domain)
                current_domain_constraints.add(variable_in_domain_constraint)
            for domain in (domains - set(domain_subset)):
                variable_notin_domain_constraint = NotInclusionConstraint(variable, domain)
                current_domain_constraints.add(variable_in_domain_constraint)
            inclusion_constraint_sets.add(ConstraintSet(current_domain_constraints))
        return inclusion_constraint_sets
