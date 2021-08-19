"""File for independent Paired groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional, cast, List, Set, FrozenSet, Sequence
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class IndependentPairedGroundings(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['VariableEquivalenceClass']]:
        """IndependentPairedGroundings is applicable if the theory is shattered
        (which it must already be to reach this rule)
        and there is a root unifying class with two variables per clause.
        Returns True and the root unifying class if applicable, and False, None otherwise."""
        unifying_classes = cnf.get_unifying_classes()
        is_root_in_cnf = lambda eq_class: eq_class.is_root_in_cnf(cnf)
        root_unifying_classes = filter(is_root_in_cnf, unifying_classes)
        for root_unifying_class in sorted(root_unifying_classes):
            if cnf.eq_class_has_two_variables(root_unifying_class):
                return True, root_unifying_class
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', root_unifying_class: 'VariableEquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentPairedGroundings and return an NNFNode"""

        # choose a clause from which to get the root variable's domain
        clause = sorted(cnf.c_clauses)[0]
        # we only need to look at one variable, because shattering means all variables in the
        # root unifying class have the same domain
        root_variable = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        domain = clause.cs.get_domain_for_variable(root_variable)

        # we substitute in two new variables for IPG
        new_x_variable = FreeVariable('X', domain)
        new_y_variable = FreeVariable('Y', domain)
        new_variables = (new_x_variable, new_y_variable)
        # we want to split the root unifying class into two,
        # each set having one variable in each clause
        root_class_a: Set['LogicalVariable'] = set()
        root_class_b: Set['LogicalVariable'] = set()
        for clause in sorted(cnf.c_clauses):
            clause_vars: FrozenSet['LogicalVariable'] = clause.bound_vars.intersection(root_unifying_class.members)
            variable_1, variable_2 = sorted(clause_vars)
            root_class_a.add(variable_1)
            root_class_b.add(variable_2)
        # we want a substitution both ways --
        # (X to root_class_a and Y to root_class_b, Y to root_class_a, X to root_class_b)
        xy_substitution_pairs_for_x = [(root_var, new_x_variable) for root_var in root_class_a]
        xy_substitution_pairs_for_y = [(root_var, new_y_variable) for root_var in root_class_b]

        yx_substitution_pairs_for_x = [(root_var, new_x_variable) for root_var in root_class_b]
        yx_substitution_pairs_for_y = [(root_var, new_y_variable) for root_var in root_class_a]

        xy_substitution = Substitution([*xy_substitution_pairs_for_x, *xy_substitution_pairs_for_y])
        yx_substitution = Substitution([*yx_substitution_pairs_for_x, *yx_substitution_pairs_for_y])

        all_x_substitution = Substitution([*xy_substitution_pairs_for_x, *yx_substitution_pairs_for_x])
        all_y_substitution = Substitution([*xy_substitution_pairs_for_y, *yx_substitution_pairs_for_y])
        new_cnf_xy = cls._substitute_root_vars(cnf, new_variables, xy_substitution)
        new_cnf_yx = cls._substitute_root_vars(cnf, new_variables, yx_substitution)
        new_cnf = new_cnf_xy.join(new_cnf_yx)

        new_x_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, all_x_substitution)
        new_y_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, all_y_substitution)
        new_xy_variables_cs = new_x_variable_cs.join(new_y_variable_cs)
        # adding the less-than constraint between the new variables
        new_variables_cs = new_xy_variables_cs.join(ConstraintSet([LessThanConstraint(new_x_variable, new_y_variable)]))
        return ForAllNode(compiler.compile(new_cnf), new_variables, new_variables_cs)

    @classmethod
    def _substitute_root_vars(cls, cnf: 'CNF', new_variables: Sequence['LogicalVariable'], substitution: 'Substitution') -> 'CNF':
        """After applying IPG we end up with a clause with free variables. This function substitutes the free variables
         for bound variables."""

        new_clauses: Set['Clause'] = set()
        # NOTE: all clauses are constrained (since must have at least two bound vars)
        for clause in cnf.c_clauses:
            new_literals = [literal.substitute(substitution) for literal in clause.literals]
            subbed_cs = clause.cs.substitute(substitution)
            if subbed_cs is None:
                raise ValueError(f"subbed_cs {subbed_cs} shouldn't be unsatisfiable")

            new_cs: 'ConstraintSet' = subbed_cs.drop_constraints_involving_only_specific_variables(new_variables)
            new_bound_vars = cast(List['LogicalVariable'],
                                  [substitution[var] for var in clause.bound_vars
                                   if not substitution[var] in new_variables])

            # add clauses back to new CNF
            if new_bound_vars == [] and new_cs == ConstraintSet([]):
                new_clauses.add(UnconstrainedClause(new_literals))
            else:
                new_clauses.add(ConstrainedClause(new_literals, new_bound_vars, new_cs))
        # shattering is preserved during this operation
        return CNF(new_clauses, shattered=cnf.shattered, subdivided=cnf.subdivided)
            
    @classmethod
    def _get_new_variable_cs(cls,
                             cnf: 'CNF',
                             root_unifying_class: 'VariableEquivalenceClass',
                             sub: 'Substitution'
                             ) -> 'ConstraintSet':
        """Return a constraint set for the new variables that has the same solutions as the root unifying variables"""
        clause = sorted(cnf.c_clauses)[0]
        # it doesn't matter which root variable we get, shattering means all solutions are the same
        root_variable = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        set_type_constraints = [sc for sc in clause.cs.set_constraints if sc.logical_term == root_variable]
        new_cs = ConstraintSet(set_type_constraints).substitute(sub)
        if new_cs is None:
            raise ValueError(f"{new_cs = } is unsatisfiable!")
        return new_cs
