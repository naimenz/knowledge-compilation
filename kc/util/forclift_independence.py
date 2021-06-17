"""This file is for translating Forclift's Independence and unify algorithms
into Python so I understand how they work and can write them more cleanly"""

from kc.data_structures import *

from typing import Optional, Tuple, List

def tryIndependentSubtheories(cnf: 'CNF'
                              ) -> Optional[Tuple['CNF', 'CNF']]:
    """They return an NNFNode but I don't have that set up yet, so I'll just
    return the independent subtheories"""

    def partition(depClauses: List['ConstrainedClause'],
                  indepClauses: List['ConstrainedClause']
                  ) -> Tuple[List['ConstrainedClause'], List['ConstrainedClause']]:
        """They use this function to construct the independent subtheories recursively"""
        if len(indepClauses) == 0:
            return depClauses, []
        # this is done with a pattern match in scala
        else:
            if len(depClauses) > 0:
                clause, rest = depClauses[0], depClauses[1:]
                indep = [other_clause for other_clause in indepClauses if clause.independent(other_clause)]
                dep = [other_clause for other_clause in indepClauses if not clause.independent(other_clause)]
                depAll, indepAll = partition(rest + dep, indep)
                return [clause] + depAll, indepAll
            else:
                return [], indepClauses

    clauses = list(cnf.clauses)
    dep, indep = partition([clauses[0]], clauses[1:])
    if len(indep) == 0:
        return None
    else:
        return CNF(dep), CNF(indep)

