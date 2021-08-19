"""File for independent single groundings compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from typing import Tuple, Optional, cast, List, Set
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class IndependentSingleGroundings(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['VariableEquivalenceClass']]:
        """IndependentSingleGroundings is applicable if the theory is shattered
        (which it must already be to reach this rule)
        and there is a root unifying class with one variable per clause.
        Returns (True, root unifying class) if applicable, and (False, None) otherwise."""
        unifying_classes = cnf.get_unifying_classes()
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
        clause = sorted(cnf.c_clauses)[0]
        # there must only be one root variable in the intersection
        root_term = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        root_variable = cast(LogicalVariable, root_term)
        domain = clause.cs.get_domain_for_variable(root_variable)

        symbol = root_variable.symbol[0]  # get just the letter, not any numbers, e.g. 'X1' -> 'X'
        new_variable = FreeVariable(symbol, domain)  

        root_substitution_pairs = sorted([(root_var, new_variable) for root_var in root_unifying_class])
        substitution = Substitution(root_substitution_pairs)
        new_variable_cs = cls._get_new_variable_cs(cnf, root_unifying_class, substitution)
        new_cnf = cls._substitute_root_vars(cnf, new_variable, substitution, new_variable_cs)
        return ForAllNode(compiler.compile(new_cnf), [new_variable], new_variable_cs)


    @classmethod
    def _substitute_root_vars(cls,
                              cnf: 'CNF',
                              new_variable: 'LogicalVariable',
                              substitution: 'Substitution',
                              new_variable_cs: 'ConstraintSet'
                              ) -> 'CNF':
        """We can't simply use the substitute method because we
        need to remove the new variable from the bound vars after substituting it.
        We also handle dropping unnecessary constraints here"""

        new_clauses: Set['Clause'] = set()
        # NOTE: (since all clauses must have at least one bound var) all clauses are constrained
        for clause in sorted(cnf.c_clauses):
            new_literals = [literal.substitute(substitution) for literal in clause.literals]
            subbed_cs = clause.cs.substitute(substitution)
            if subbed_cs is None:
                raise ValueError(f"subbed_cs {subbed_cs} shouldn't be unsatisfiable")
            # because many constraints are now taken care of in the new forall node, we need to remove them
            # from the clause. e.g. if we had \forall X, Y, X != Y: p(X, Y) and ISG removed X,
            # then we would want to remove the X != Y constraint from the new clause, giving
            # (\forall X, X != Y) [\forall Y: p(X, Y)]
            redundant_set_constraints = set(c for c in subbed_cs.set_constraints if c.logical_term == new_variable)
            new_cs: 'ConstraintSet' = ConstraintSet(subbed_cs.constraints
                                                    - new_variable_cs.constraints
                                                    - redundant_set_constraints)
            new_bound_vars = cast(List['LogicalVariable'],
                                  [substitution[var] for var in clause.bound_vars if substitution[var] != new_variable])

            # add clauses back to new CNF
            if new_bound_vars == [] and new_cs == ConstraintSet([]):  # (\forall X, X != Y) [p(X, Y)] -> add p(X,Y)
                new_clauses.add(UnconstrainedClause(new_literals))
            else:  # (\forall X, X != Y) [\forall Y: p(X, Y)] -> add [\forall Y: p(X, Y)]
                new_clauses.add(ConstrainedClause(new_literals, new_bound_vars, new_cs))
        # shattering is preserved during this operation
        return CNF(new_clauses, shattered=cnf.shattered, subdivided=cnf.subdivided, names=None)
            
    @classmethod
    def _get_new_variable_cs(cls,
                             cnf: 'CNF',
                             root_unifying_class: 'VariableEquivalenceClass',
                             sub: 'Substitution',
                             domain: 'ProperDomain'
                             ) -> 'ConstraintSet':
        """Return a constraint set for the new variable that has the same solutions as the root unifying variables"""
        # we only loop once to get a clause from the cnf
        clause = sorted(cnf.c_clauses)[0]
        # must only be one root variable in the intersection
        root_term = sorted(root_unifying_class.members.intersection(clause.bound_vars))[0]
        root_variable = cast(LogicalVariable, root_term)
        # for checking if the logical constraints are redundant
        domain = clause.cs.get_domain_for_variable(root_variable)
        # NOTE : to get the WFOMI computation right, we want to include the inequality with the BIGGER domain by
        # looking at the domain of the other variable, e.g. \forall X \in D, \forall \y in D_top, X!=Y
        # then D_top being a subdomain we include the inequality X!=Y in D
        # since we are in the two variable fragment, there are 0 or 1 remaining variables after ISG
        other_variable_set = clause.bound_vars - {root_variable}
        if len(other_variable_set) > 0:
            next_variable = get_element_of_set(other_variable_set)
            next_domain = clause.cs.get_domain_for_variable(next_variable)
            next_domain_is_subset = next_domain.is_strict_subset_of(domain)
        else:
            # there is no next domain, so in some sense it is a subset, but really this is just for convenience
            next_domain_is_subset = True
        set_constraints = [sc for sc in clause.cs.set_constraints if sc.logical_term == root_variable]

        # only include logical constraints if this is the larger domain
        if next_domain_is_subset:
            logical_constraints = [lc for lc in clause.cs.logical_constraints
                                   if cls._is_valid_logical_constraint(lc, root_variable, domain)]
        else:
            logical_constraints = []
        new_cs = ConstraintSet([*set_constraints, *logical_constraints]).substitute(sub)
        if new_cs is None:
            raise ValueError(f"{new_cs = } is unsatisfiable!")
        return new_cs

    @classmethod
    def _is_valid_logical_constraint(cls,
                                     constraint: 'LogicalConstraint',
                                     root_variable: 'LogicalVariable',
                                     domain: 'ProperDomain') -> bool:
        left_term = constraint.left_term
        right_term = constraint.right_term
        if left_term == root_variable:
            if isinstance(right_term, FreeVariable) and right_term.domain.intersect_with(domain) != EmptyDomain():
                return True
        elif right_term == root_variable:
            if isinstance(left_term, FreeVariable) and left_term.domain.intersect_with(domain) != EmptyDomain():
                return True
        return False

