"""
This file contains various useful functions for working with the formulas defined in data_structures
"""

from kc.data_structures import *
from kc.util import *

from copy import deepcopy
from itertools import chain, product

from typing import List, Tuple, Set, Dict, cast, FrozenSet, Optional
from typing import TypeVar

VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']
C = TypeVar('C', bound='ConstrainedClause') # for functions that can operate on more than one type of clause

def get_constrained_atoms(clause: 'ConstrainedClause') -> List['ConstrainedAtom']:
    """For a given clause, return a list of all the constrained atoms in the clause"""
    constrained_atoms = []
    for literal in clause.unconstrained_clause.literals:
        constrained_atom = build_constrained_atom_from_literal(literal, clause)
        constrained_atoms.append(constrained_atom)
    return constrained_atoms


def build_constrained_atom_from_literal(literal: 'Literal', clause: 'ConstrainedClause') -> 'ConstrainedAtom':
    """Build a constrained atom given a specific literal and its parent clause"""
    atom = literal.atom 
    positive_literal = Literal(atom, True)
    # constrained atoms subclass ConstrainedClause, so need to be built of an
    # unconstrained clause
    unconstrained_atom = UnconstrainedClause([positive_literal])
    constrained_atom = ConstrainedAtom(unconstrained_atom, clause.bound_vars, clause.cs)
    return constrained_atom


def get_constrained_literals(clause: 'ConstrainedClause') -> List['UnitClause']:
    """For a given clause, return a list of all the constrained literals in the clause"""
    constrained_literals = []
    for literal in clause.unconstrained_clause.literals:
        constrained_literal = build_constrained_literal_from_literal(literal, clause)
        constrained_literals.append(constrained_literal)
    return constrained_literals


def build_constrained_literal_from_literal(literal: 'Literal', clause: 'ConstrainedClause') -> 'UnitClause':
    """Build a constrained literal given a specific literal and its parent clause"""
    unconstrained_literal = UnconstrainedClause([literal])
    constrained_literal = UnitClause(unconstrained_literal, clause.bound_vars, clause.cs)
    return constrained_literal
