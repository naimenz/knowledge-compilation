"""File for independence compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import *

from typing import Tuple, Optional, Sequence
StoredClauses = Tuple[Sequence['Clause'], Sequence['Clause']]

class Independence(KCRule):
    @staticmethod
    def is_applicable(delta: 'CNF') -> Tuple[bool, Optional[StoredClauses]]:
        """Independence is applicable if the theory can be divided into
        two subtheories such that the subtheories make up the whole theory and are
        independent.
        NOTE: This will work with domai nvariables when clauses_independent does"""
        clauses = list(delta.clauses)
        potential_subtheory, other_potential_subtheory = Independence._partition([clauses[0]], clauses[1:])
        # if all clauses have been moved into potential_subtheory, then we are back where we started!
        if len(other_potential_subtheory) == 0:
            return False, None
        else:
            return True, (potential_subtheory, other_potential_subtheory)

    @staticmethod
    def apply(delta: 'CNF', stored_data: StoredClauses) -> 'NNFNode':
        """Apply Independence and return an NNFNode"""
        raise NotImplementedError('Independence.apply not implemented')


    @staticmethod
    def _partition(potential_subtheory: List['ConstrainedClause'],
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
                clause, rest = potential_subtheory[0], other_clauses[1:]
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
                potential_subtheory, other_clauses = Independence._partition(rest + new_potential_subtheory, new_other_clauses)
                return [clause] + potential_subtheory, other_clauses
            else:
                return [], other_clauses
