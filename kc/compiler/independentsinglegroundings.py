"""File for independent single groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from typing import Tuple, Optional, cast, List, Set
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
        for root_unifying_class in sorted(root_unifying_classes):
            if cnf.eq_class_has_one_variable(root_unifying_class):
                return True, root_unifying_class
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', root_unifying_class: 'VariableEquivalenceClass', compiler: 'Compiler') -> 'NNFNode':
        """Apply IndependentSingleGroundings and return an NNFNode"""
        
        # choose a clause to get the root variable from
        for clause in sorted(cnf.c_clauses): break
        # must only be one root variable in the intersection
        root_term = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        root_variable = cast(LogicalVariable, root_term)  # hack for type checking
        domain = clause.cs.get_domain_for_variable(root_variable)

        symbol = root_variable.symbol[0]  # get just the letter, not any numbers
        new_variable = FreeVariable(symbol, domain)  

        root_substitution_pairs = sorted([(root_var, new_variable) for root_var in root_unifying_class])
        substitution = Substitution(root_substitution_pairs)
        new_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, substitution)
        new_cnf = cls._substitute_root_vars(cnf, new_variable, substitution, new_variable_cs)
        return ForAllNode(compiler.compile(new_cnf), [new_variable], new_variable_cs)


    @classmethod
    def _substitute_root_vars(cls, cnf: 'CNF', new_variable: 'LogicalVariable', sub: 'Substitution', new_variable_cs: 'ConstraintSet') -> 'CNF':
        """We can't simply use the substitute method because we
        need to remove the new variable from the bound vars after substituting it.
        We also handle dropping unnecessary constraints here"""

        new_clauses: Set['Clause'] = set()
        # NOTE: all clauses are constrained (since must have at least one bound var)
        for clause in sorted(cnf.c_clauses):
            new_literals = [literal.substitute(sub) for literal in clause.literals]
            subbed_cs = clause.cs.substitute(sub)
            if subbed_cs is None:
                raise ValueError(f"subbed_cs {subbed_cs} shouldn't be unsatisfiable")
            # NOTE TODO: Trying out removing the specific constraints used rather than all with the variable
            # and we also remove any set constraints that involve the variable 
            redundant_set_constraints = set(c for c in subbed_cs.set_constraints if c.logical_term == new_variable)
            new_cs: 'ConstraintSet' = ConstraintSet(subbed_cs.constraints - new_variable_cs.constraints - redundant_set_constraints)
            _new_bound_vars = [sub[var] for var in clause.bound_vars if sub[var] != new_variable]
            new_bound_vars = cast(List['LogicalVariable'], _new_bound_vars) # hack for type checking

            # handle making unconstrained clauses
            if new_bound_vars == [] and new_cs == ConstraintSet([]):
                new_clauses.add(UnconstrainedClause(new_literals))
            else:
                new_clauses.add(ConstrainedClause(new_literals, new_bound_vars, new_cs))
        # shattering is preserved during this operation
        # DEBUG shouldnt need to rename here
        return CNF(new_clauses, shattered = cnf.shattered, names=None)
            
    @classmethod
    def _get_new_variable_cs(cls,
                             cnf: 'CNF',
                             root_unifying_class: 'VariableEquivalenceClass',
                             sub: 'Substitution'
                             ) -> 'ConstraintSet':
        """Return a constraint set for the new variable that has the same solutions as the root unifying variables"""
        # we only loop once to get a clause from the cnf
        for clause in sorted(cnf.c_clauses): break
        # must only be one root variable in the intersection
        root_term = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        root_variable = cast(LogicalVariable, root_term)  # hack for type checking
        # for checking if the logical constraints are redundant
        domain = clause.cs.get_domain_for_variable(root_variable) 
        # NOTE TODO DEBUG: trying out including the inequality with the BIGGER domain by looking at the domain of the other var
        # since we are in the two variable fragment, there are 0 or 1 remaining variables after ISG
        other_variable_set = clause.bound_vars - set([root_variable])
        if len(other_variable_set) > 0:
            next_variable = get_element_of_set(other_variable_set)
            next_domain = clause.cs.get_domain_for_variable(next_variable)
            next_domain_is_subset = next_domain.is_strict_subset_of(domain)
        else:
            # there is no next domain, so in some sense it is a subset, but really this is just for convenience
            next_domain_is_subset = True
        set_constraints = [sc for sc in clause.cs.set_constraints if sc.logical_term == root_variable]

        # NOTE DEBUG TODO: Trying a small function to check whether we should include this constraint yet
        def is_valid_logical_constraint(constraint: 'LogicalConstraint') -> bool:
            lt = constraint.left_term
            rt = constraint.right_term
            # only include logical constraints if this is the larger domain
            if next_domain_is_subset:
                if lt == root_variable:
                    if isinstance(rt, FreeVariable) and rt.domain.intersect_with(domain) != EmptyDomain():
                        return True
                elif rt == root_variable:
                    if isinstance(lt, FreeVariable) and lt.domain.intersect_with(domain) != EmptyDomain():
                        return True
            return False

        logical_constraints = [lc for lc in clause.cs.logical_constraints if is_valid_logical_constraint(lc)]
        new_cs = ConstraintSet([*set_constraints, *logical_constraints]).substitute(sub)
        if new_cs is None:
            raise ValueError(f"{new_cs = } is unsatisfiable!")
        return new_cs


        
