"""File for independence compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import *

from typing import Tuple, Optional, Sequence, List, Any
StoredCNFs = Tuple['CNF', 'CNF']

class Independence(KCRule):
    @classmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional[StoredCNFs]]:
        """Independence is applicable if the theory can be divided into
        two subtheories such that the subtheories make up the whole theory and are
        independent.
        NOTE: This will work with domain variables when clauses_independent does"""
        # TODO: can partition be rewritten to take sets?
        clauses = list(delta.clauses)
        subtheory, other_subtheory = cls._partition([clauses[0]], clauses[1:])
        # if all clauses have been moved into subtheory, then we are back where we started!
        if len(other_subtheory) == 0:
            return False, None
        else:
            return True, (CNF(subtheory), CNF(other_subtheory))

    @classmethod
    # NOTE: We annotate compiler as 'Any' to avoid circularly importing the Compiler
    def apply(cls, delta: 'CNF', sub_cnfs: StoredCNFs, compiler: Any) -> 'NNFNode':
        """Apply Independence and return an NNFNode (in this case an AndNode)"""
        return AndNode(compiler.compile(sub_cnfs[0]), compiler.compile(sub_cnfs[1]))

    @classmethod
    def _partition(cls, potential_subtheory: List['ConstrainedClause'],
                  other_clauses: List['ConstrainedClause']
                  ) -> Tuple[List['ConstrainedClause'], List['ConstrainedClause']]:
        """They use this function to construct the independent subtheories recursively.
        We start with two guesses for independent subtheories.
        We then iteratively move all the clauses that are dependent with the potential_subtheory 
        over to it, and keep the rest separately.
        """
        if len(other_clauses) == 0:
            return potential_subtheory, []
        # this is done with a pattern match in scala
        else:
            if len(potential_subtheory) > 0:
                # TODO:  check this line with Forclift
                clause, rest = potential_subtheory[0], potential_subtheory[1:]
                new_potential_subtheory = []
                new_other_clauses = []

                # separate out the clauses into whether they are dependent with the potential_subtheory or not
                for other_clause in other_clauses:
                    independent = clauses_independent(clause, other_clause)
                    if not independent:
                        new_potential_subtheory.append(other_clause)
                    else:
                        new_other_clauses.append(other_clause)

                # recurse
                potential_subtheory, other_clauses = cls._partition(rest + new_potential_subtheory, new_other_clauses)
                return [clause] + potential_subtheory, other_clauses
            else:
                return [], other_clauses
