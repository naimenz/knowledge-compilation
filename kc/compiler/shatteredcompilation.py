"""File for shatter compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import powerset, partition_set

from itertools import product
from functools import reduce

from typing import Tuple, Set, Iterable, Sequence, FrozenSet
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
        shattered_clauses_list = [cls.shatter_clause(clause, terms, domains) for clause in cnf.clauses]
        empty_set: Set['ConstrainedClause'] = set() # hack for type checking
        flattened_shattered_clauses = empty_set.union(*shattered_clauses_list)
        propagated_shattered_clauses = [c.propagate_equality_constraints() for c in flattened_shattered_clauses]
        return compiler.compile(CNF(propagated_shattered_clauses, shattered=True, subdivided=cnf.subdivided))

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
        # NOTE TODO: I'm assuming here that we only look at bound variables, although this is not quite what the PhD says
        # literal_variables = clause.literal_variables
        literal_variables = clause.bound_vars
        # CS_A
        shatter_var_constraints = cls._build_shatter_var_constraints(literal_variables, terms, domains)
        # CS_B
        literal_inequality_constraints = cls._build_literal_inequality_constraints(clause)
        # CS
        all_final_constraints = product([clause.cs], shatter_var_constraints, literal_inequality_constraints)
        final_constraints = map(cls._merge_constraint_sets, all_final_constraints)

        # only return the satisfiable ones
        satisfiable_constraints = [cs for cs in final_constraints if cs.is_satisfiable()]
        return set(ConstrainedClause(clause.literals, clause.bound_vars, cs) for cs in satisfiable_constraints)

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
        final_constraints = set(cs_eq.join(cs_in) for cs_eq, cs_in in product_of_cs)
        # only return satisfiable constraint sets
        return set(cs for cs in final_constraints if cs.is_satisfiable())

    @classmethod
    def _build_literal_inequality_constraints(cls, clause: 'ConstrainedClause') -> Set['ConstraintSet']:
        """Build all of the equality and inequality constraints between variables in each literal
        (CS_B from the PhD)."""
        individual_literal_inequality_constraints = []
        for literal in sorted(clause.literals):
            bound_literal_variables = set(clause.bound_vars.intersection(literal.variables))
            literal_constraint_sets = cls._build_literal_constraint_sets(bound_literal_variables)
            individual_literal_inequality_constraints.append(literal_constraint_sets)

        all_literal_inequality_constraints = product(*individual_literal_inequality_constraints)
        # merge each combination into a single constraint set 
        literal_inequality_constraints = map(cls._merge_constraint_sets, all_literal_inequality_constraints)
        return set(literal_inequality_constraints)

    @classmethod
    def _build_literal_constraint_sets(cls, bound_literal_variables: Set['LogicalVariable']) -> Set['ConstraintSet']:
        """Build the (in)equality between variable constraints for this literal 
        (CS_B^i from the PhD)"""
        literal_constraint_sets: Set['ConstraintSet'] = set()
        for partition in partition_set(bound_literal_variables):
            partition_constraint_set = cls._build_partition_constraint_set(partition)
            literal_constraint_sets.add(partition_constraint_set)
        return literal_constraint_sets

    @classmethod
    def _build_partition_constraint_set(cls, partition: Set[FrozenSet['LogicalVariable']]) -> 'ConstraintSet':
        """Build the constraint set for a particular partition of the logical variables"""
        equality_constraints: Set['EqualityConstraint'] = set()
        inequality_constraints: Set['InequalityConstraint'] = set()
        for subset in sorted(partition): 
            # there is some duplicated effort here, but it is still correct because of sets
            # TODO: for performance, reduce looping
            for var in sorted(subset):
                for other_var in sorted(subset - set([var])):
                    equality_constraints.add(EqualityConstraint(var, other_var))

                for other_subset in sorted(partition - set([subset])):
                    for other_var in sorted(other_subset):
                        inequality_constraints.add(InequalityConstraint(var, other_var))
        # we then combine all of these into a single constraint
        partition_constraint_set = ConstraintSet([*equality_constraints, *inequality_constraints])
        return partition_constraint_set

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
        """Merge a collection of constraint sets into one"""
        return reduce(lambda cs1, cs2: cs1.join(cs2), cs_iterable, ConstraintSet([]))

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
        for free_variable in sorted(free_variables):
            var_eq_constraint = EqualityConstraint(variable, free_variable)
            equality_constraint_sets.add(ConstraintSet([var_eq_constraint]))
            not_equal_constraints.add(~var_eq_constraint)
        for constant in sorted(constants):
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
        for domain_subset in sorted(all_domain_subsets):
            current_domain_constraints: Set['SetConstraint'] = set()
            for domain in sorted(domain_subset):
                variable_in_domain_constraint = InclusionConstraint(variable, domain)
                current_domain_constraints.add(variable_in_domain_constraint)
            for domain in sorted(domains - set(domain_subset)):
                variable_notin_domain_constraint = NotInclusionConstraint(variable, domain)
                current_domain_constraints.add(variable_notin_domain_constraint)
            inclusion_constraint_sets.add(ConstraintSet(current_domain_constraints))
        return inclusion_constraint_sets
