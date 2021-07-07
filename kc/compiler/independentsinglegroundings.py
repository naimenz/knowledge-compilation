"""File for independent single groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional, cast, List
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class IndependentSingleGroundings(KCRule):
    
    # TODO: Try checking ISG and IPG together
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['VariableEquivalenceClass']]:
        """IndependentSingleGroundings is applicable if the theory is shattered
        (which it must already be to reach this rule)
        and there is a root unifying class with one variable per clause.
        Returns True and the root unifying class if applicable, and False, None otherwise."""
        unifying_classes = cnf.get_unifying_classes()
        # writing a little function to check if root in this particular cnf
        # TODO: Think about making this more readable or faster
        is_root_in_cnf = lambda eq_class: eq_class.is_root_in_cnf(cnf)
        root_unifying_classes = filter(is_root_in_cnf, unifying_classes)
        for root_unifying_class in root_unifying_classes:
            if cnf.eq_class_has_one_variable(root_unifying_class):
                return True, root_unifying_class
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', root_unifying_class: 'VariableEquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentSingleGroundings and return an NNFNode"""
        new_variable = cnf.get_new_logical_variable('X')
        root_substitution_pairs = [(root_var, new_variable) for root_var in root_unifying_class]
        substitution = Substitution(root_substitution_pairs)
        new_cnf = cls._substitute_root_vars(cnf, new_variable, substitution)
        new_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, substitution)
        return ForAllNode(compiler.compile(new_cnf), [new_variable], new_variable_cs)


    @classmethod
    def _substitute_root_vars(cls, cnf: 'CNF', new_variable: 'LogicalVariable', sub: 'Substitution') -> 'CNF':
        """We can't simply use the substitute method because we
        need to remove the new variable from the bound vars after substituting it """

        new_clauses = set()
        # NOTE: all clauses are constrained (since must have at least one bound var)
        for clause in cnf.c_clauses:
            new_literals = [literal.substitute(sub) for literal in clause.literals]
            new_cs = clause.cs.substitute(sub).drop_constraints_involving_only_specific_variables([new_variable])
            _new_bound_vars = [sub[var] for var in clause.bound_vars if sub[var] != new_variable]
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
        """Return a constraint set for the new variable that has the same solutions as the root unifying variables"""
        # we only loop once to get a clause from the cnf
        for clause in cnf.c_clauses: break
        # must only be one root variable in the intersection
        root_variable = list(root_unifying_class.members.intersection(clause.bound_vars))[0]
        set_constraints = [sc for sc in clause.cs.set_constraints if sc.logical_term == root_variable]
        new_cs = ConstraintSet(set_constraints)
        return new_cs.substitute(sub)


        
