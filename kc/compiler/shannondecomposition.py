"""File for shannon decomposition compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from typing import Tuple, Optional, List
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler
class ShannonDecomposition(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional['ConstrainedAtom']]:
        """ShannonDecomposition is applicable if the theory contains
        an atom without bound logical variables.
        Returns True plus the atom if applicable, and False, None otherwise"""
        for clause in sorted(cnf.clauses):
            for c_atom in sorted(clause.get_constrained_atoms()):
                overlap = c_atom.bound_vars.intersection(set(c_atom.atom.terms)) 
                if len(overlap) == 0:
                    return True, c_atom
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', unbound_atom: 'ConstrainedAtom', compiler: 'Compiler') -> 'NNFNode':
        """Apply ShannonDecomposition and return an NNFNode"""
        # NOTE: these 'literals' are technically a single-item set of literals
        true_literal = get_element_of_set(unbound_atom.literals)
        true_branch = cnf.join(CNF([UnconstrainedClause([true_literal])]))
        if not true_literal.is_smt():
            false_literal = ~true_literal
            false_branch = cnf.join(CNF([UnconstrainedClause([false_literal])]))
        # with SMT atoms, we don't use negated literals. Instead, we take the outer ranges 
        # (e.g. ¬[10 < age(X) < 20] goes to -inf < age(X) < 10 v 20 < age(X) < inf)
        else:
            old_predicate = true_literal.atom.predicate
            new_low_predicate = SMTPredicate(old_predicate.name, old_predicate.arity, float('-inf'), old_predicate.lower_bound)
            new_high_predicate = SMTPredicate(old_predicate.name, old_predicate.arity, old_predicate.upper_bound, float('-inf'))
            terms = true_literal.atom.terms
            false_literals = [Literal(Atom(new_low_predicate, terms)), Literal(Atom(new_high_predicate, terms))]
            false_branch = cnf.join(CNF([UnconstrainedClause(false_literals)]))
        # manually setting flags because cnf.join cannot
        true_branch.shattered = cnf.shattered
        true_branch.subdivided = cnf.subdivided
        false_branch.shattered = cnf.shattered
        false_branch.subdivided = False  # since we have new predicates, we may not be subdivided

        return OrNode(compiler.compile(true_branch), compiler.compile(false_branch))

