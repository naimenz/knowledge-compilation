"""
This is a file for parsing input to knowledge compilation.
"""

from kc.data_structures import *

from typing import List

def make_auxiliary_predicate_for_clause(clause: 'ConstrainedClause', auxiliary_name: str='a1') -> List['ConstrainedClause']:
    """Given a CNF clause, construct an auxiliary predicate for the clause and return a CNF
    that is equivalent containing that auxiliary predicate

    To enforce that the aux pred (aux) is equivalent to the clause (clause), we need (aux => clause) and (clause => aux).
    Equivalently, we need (¬clause => ¬aux)  and (clause => aux)
    We'll call (¬clause => ¬aux) the first branch. It is straightforward and gives (clause v ¬aux).
    The second branch gives a conjunction of a number of formulas of the form (aux v ¬l) for each literal l in clause """
    # NOTE: I'm assuming the aux pred takes as many arguments as the clause has bound variables
    aux_terms = list(clause.bound_vars)
    auxiliary_predicate = Predicate(auxiliary_name, len(aux_terms))
    aux_literal = Literal(Atom(auxiliary_predicate, aux_terms))

    first_branch = [ConstrainedClause(list(clause.literals) + [~aux_literal], clause.bound_vars, clause.cs)]

    second_branch = [ConstrainedClause([literal, aux_literal], clause.bound_vars, clause.cs) for literal in clause.literals]

    return first_branch + second_branch

