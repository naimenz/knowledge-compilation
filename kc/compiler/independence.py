"""File for independence compilation rule"""

from kc.data_structures import *
from kc.compiler import KCRule

from typing import Tuple, Optional, List
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class Independence(KCRule):
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional[Tuple['CNF', 'CNF']]]:
        """Independence is applicable if the theory can be divided into two subtheories
         such that the subtheories make up the whole theory and are independent."""
        if len(cnf.clauses) == 1:
            return False, None  # need at least 2 clauses to make non-empty independent sets
        cnf_clauses = sorted(cnf.clauses)
        subtheory, other_subtheory = cls._partition([cnf_clauses[0]], cnf_clauses[1:])
        # if all clauses have been moved into subtheory by cls._partition, no theories are independent!
        if len(other_subtheory) == 0:
            return False, None
        else:
            # if the parent was shattered, so are the children
            return True, (CNF(subtheory, shattered=cnf.shattered, subdivided=cnf.subdivided),
                          CNF(other_subtheory, shattered=cnf.shattered, subdivided=cnf.subdivided))

    @classmethod
    def apply(cls, cnf: 'CNF', sub_cnfs: Tuple['CNF', 'CNF'], compiler: 'Compiler') -> 'NNFNode':
        """Apply Independence and return an NNFNode (in this case an AndNode)"""
        return AndNode(compiler.compile(sub_cnfs[0]), compiler.compile(sub_cnfs[1]))

    @classmethod
    def _partition(cls, potential_subtheory: List['Clause'],
                   other_clauses: List['Clause']
                   ) -> Tuple[List['Clause'], List['Clause']]:
        """Used to construct the independent subtheories recursively.
        We start with two guesses for independent subtheories.
        We then iteratively move all the clauses that are dependent with the potential_subtheory 
        over into one list, and keep the rest separately.
        """
        if len(other_clauses) == 0:
            return potential_subtheory, []
        else:
            if len(potential_subtheory) > 0:
                clause, rest = potential_subtheory[0], potential_subtheory[1:]
                new_potential_subtheory = []
                new_other_clauses = []

                # separate out the clauses into whether they are dependent with the potential_subtheory or not
                for other_clause in other_clauses:
                    independent = clause.is_independent_from_other_clause(other_clause)
                    if not independent:
                        new_potential_subtheory.append(other_clause)
                    else:
                        new_other_clauses.append(other_clause)

                # recurse
                potential_subtheory, other_clauses = cls._partition(rest + new_potential_subtheory, new_other_clauses)
                return [clause] + potential_subtheory, other_clauses
            else:
                return [], other_clauses
