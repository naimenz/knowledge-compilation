"""
This file contains various useful functions for working with the formulas defined in data_structures
"""

from kc.data_structures.cnf import *
from kc.data_structures.literals import *
from kc.data_structures.constraints import *

from typing import List


# def get_predicates(cnf: 'FO_CNF') -> list['Predicate']: 
#     pred = Predicate('fun', 1)
#     return [pred]

def get_constrained_atoms(clause: 'ConstrainedClause') -> List['ConstrainedAtom']:
    """For a given clause, return a list of all the constrained atoms in the clause"""
    constrained_atoms = []
    for literal in clause.unconstrained_clause.literals:
        atom = literal.atom 
        positive_literal = Literal(atom, True)
        # constrained atoms subclass ConstrainedClause, so need to be built of an
        # unconstrained clause
        unconstrained_atom = UnconstrainedClause([positive_literal])
        constrained_atom = ConstrainedAtom(unconstrained_atom, clause.bound_vars, clause.cs)
        constrained_atoms.append(constrained_atom)
    return constrained_atoms


if __name__ == '__main__':
    preds = [Predicate('smokes', 2), Predicate('friends', 2), Predicate('fun', 1)]
    constants = [Constant('a'), Constant('b'), Constant('c')]
    variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]

    atoms = [Atom(preds[0], [variables[0], variables[1]]),
             Atom(preds[1], [variables[1], constants[0]]),
             Atom(preds[2], [constants[0]])
            ]

    literals = [Literal(atoms[0], True), Literal(atoms[1], False), Literal(atoms[2], False)]

    constraint_set = ConstraintSet([EqualityConstraint(variables[0], constants[0]),
                      InequalityConstraint(variables[1], constants[1])
                      ])

    u_clause = UnconstrainedClause(literals)
    c_clause = ConstrainedClause(u_clause, variables[:2], constraint_set)

    c_atoms = get_constrained_atoms(c_clause)
    for c_atom in c_atoms:
        print(c_atom)

