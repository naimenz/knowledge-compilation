"""File for independent Paired groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from typing import Tuple, Optional, cast, List, Iterable, Set, FrozenSet, Sequence
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
        # writing a little function to check if root in this particular cnf
        is_root_in_cnf = lambda eq_class: eq_class.is_root_in_cnf(cnf)
        root_unifying_classes = filter(is_root_in_cnf, unifying_classes)
        for root_unifying_class in root_unifying_classes:
            if cnf.eq_class_has_two_variables(root_unifying_class):
                return True, root_unifying_class
        return False, None


    @classmethod
    def apply(cls, cnf: 'CNF', root_unifying_class: 'VariableEquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentPairedGroundings and return an NNFNode"""
        # we substitute in two new variables for Paired
        new_Y_variable = cnf.get_new_logical_variable('Y')
        new_Z_variable = cnf.get_new_logical_variable('Z')
        new_variables = (new_Y_variable, new_Z_variable)
        # we want to split the root unifying class into two,
        # each set having one variable in each clause
        root_class_a: Set['LogicalVariable'] = set()
        root_class_b: Set['LogicalVariable'] = set()
        # TODO: compare speed with list comps
        for clause in cnf.c_clauses: 
            clause_vars: FrozenSet['LogicalVariable'] = clause.bound_vars.intersection(root_unifying_class.members)
            variable_1, variable_2 = tuple(clause_vars)
            root_class_a.add(variable_1)
            root_class_b.add(variable_2)
        # we want a substitution both ways -- (Y to U_a and Z to U_b, Z to U-a, Y to U_b in the PhD)
        YZ_substitution_pairs_for_Y = [(root_var, new_Y_variable) for root_var in root_class_a]
        YZ_substitution_pairs_for_Z = [(root_var, new_Z_variable) for root_var in root_class_b]

        ZY_substitution_pairs_for_Y = [(root_var, new_Y_variable) for root_var in root_class_b]
        ZY_substitution_pairs_for_Z = [(root_var, new_Z_variable) for root_var in root_class_a]

        YZ_substitution = Substitution([*YZ_substitution_pairs_for_Y, *YZ_substitution_pairs_for_Z])
        ZY_substitution = Substitution([*ZY_substitution_pairs_for_Y, *ZY_substitution_pairs_for_Z])

        all_Y_substitution = Substitution([*YZ_substitution_pairs_for_Y, *ZY_substitution_pairs_for_Y])
        all_Z_substitution = Substitution([*YZ_substitution_pairs_for_Z, *ZY_substitution_pairs_for_Z])
        new_cnf_YZ = cls._substitute_root_vars(cnf, new_variables, YZ_substitution)
        new_cnf_ZY = cls._substitute_root_vars(cnf, new_variables, ZY_substitution)
        new_cnf = new_cnf_YZ.join(new_cnf_ZY)

        new_Y_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, all_Y_substitution)
        new_Z_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, all_Z_substitution)
        new_YZ_variables_cs = new_Y_variable_cs.join(new_Z_variable_cs)
        # adding the less-than constraint between the new variables
        new_variables_cs = new_YZ_variables_cs.join(ConstraintSet([LessThanConstraint(new_Y_variable, new_Z_variable)]))
        return ForAllNode(compiler.compile(new_cnf), new_variables, new_variables_cs)

    @classmethod
    def _substitute_root_vars(cls, cnf: 'CNF', new_variables: Sequence['LogicalVariable'], sub: 'Substitution') -> 'CNF':
        """We can't simply use the apply_substitution method because we
        need to remove the new variables from the bound vars after substituting them """

        new_clauses = set()
        # NOTE: all clauses are constrained (since must have at least two bound vars)
        for clause in cnf.c_clauses:
            new_literals = [literal.substitute(sub) for literal in clause.literals]
            new_cs = clause.cs.substitute(sub).drop_constraints_involving_only_these_variables(new_variables)
            _new_bound_vars = [sub[var] for var in clause.bound_vars if not sub[var] in new_variables]
            new_bound_vars = cast(List['LogicalVariable'], _new_bound_vars) # hack for type checking

            new_clauses.add(ConstrainedClause(new_literals, new_bound_vars, new_cs))
        # shattering is preserved during this operation
        return CNF(new_clauses, shattered = cnf.shattered)
            
    @classmethod
    def _get_new_variable_cs(cls,
                             cnf: 'CNF',
                             root_unifying_class: 'VariableEquivalenceClass',
                             sub: 'Substitution'
                             ) -> 'ConstraintSet':
        """Return a constraint set for the new variables that has the same solutions as the root unifying variables"""
        # we only loop once to get a clause from the cnf
        for clause in cnf.c_clauses: break
        # it doesn't matter which root variable we get, shattering means all sols are the same
        root_variable = get_element_of_set(root_unifying_class.members.intersection(clause.bound_vars))
        set_constraints = [sc for sc in clause.cs.set_constraints if sc.logical_term == root_variable]
        new_cs = ConstraintSet(set_constraints)
        return new_cs.substitute(sub)


        
