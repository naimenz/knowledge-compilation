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
        for root_unifying_class in sorted(root_unifying_classes):
            if cnf.eq_class_has_two_variables(root_unifying_class):
                return True, root_unifying_class
        return False, None


    @classmethod
    def apply(cls, cnf: 'CNF', root_unifying_class: 'VariableEquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentPairedGroundings and return an NNFNode"""


        # choose a clause to get the root variable's domain from
        for clause in sorted(cnf.c_clauses): break
        # must only be one root variable in the intersection
        root_term = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        root_variable = cast(LogicalVariable, root_term)  # hack for type checking
        symbol = root_variable.symbol[0]  # get just the letter, not any numbers
        domain = clause.cs.get_domain_for_variable(root_variable)

        # we substitute in two new variables for Paired
        # NOTE: The variable names contain Y and Z in analogy with the PhD,
        # but we actually want our variables to be named X and Y
        new_Y_variable = FreeVariable('X', domain)  
        new_Z_variable = FreeVariable('Y', domain)
        new_variables = (new_Y_variable, new_Z_variable)
        # we want to split the root unifying class into two,
        # each set having one variable in each clause
        root_class_a: Set['LogicalVariable'] = set()
        root_class_b: Set['LogicalVariable'] = set()
        # TODO: compare speed with list comps
        for clause in sorted(cnf.c_clauses):
            clause_vars: FrozenSet['LogicalVariable'] = clause.bound_vars.intersection(root_unifying_class.members)
            variable_1, variable_2 = sorted(clause_vars)
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

        new_clauses: Set['Clause'] = set()
        # NOTE: all clauses are constrained (since must have at least two bound vars)
        for clause in cnf.c_clauses:
            new_literals = [literal.substitute(sub) for literal in clause.literals]
            subbed_cs = clause.cs.substitute(sub)
            if subbed_cs is None:
                raise ValueError(f"subbed_cs {subbed_cs} shouldn't be unsatisfiable")
            new_cs: 'ConstraintSet' = subbed_cs.drop_constraints_involving_only_specific_variables(new_variables)
            _new_bound_vars = [sub[var] for var in clause.bound_vars if not sub[var] in new_variables]
            new_bound_vars = cast(List['LogicalVariable'], _new_bound_vars) # hack for type checking

            # handle making unconstrained clauses
            if new_bound_vars == [] and new_cs == ConstraintSet([]):
                new_clauses.add(UnconstrainedClause(new_literals))
            else:
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
        for clause in sorted(cnf.c_clauses): break
        # it doesn't matter which root variable we get, shattering means all sols are the same
        root_variable = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        set_constraints = [sc for sc in clause.cs.set_constraints if sc.logical_term == root_variable]
        new_cs = ConstraintSet(set_constraints).substitute(sub)
        if new_cs is None:
            raise ValueError(f"{new_cs = } is unsatisfiable!")
        return new_cs


        
