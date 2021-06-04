"""
This file contains various useful functions for working with the formulas defined in data_structures
"""

from kc.data_structures.cnf import *
from kc.data_structures.literals import *
from kc.data_structures.constraints import *

from copy import deepcopy
from itertools import chain

from typing import List, Tuple, Set, Dict


def get_constrained_atoms(clause: 'ConstrainedClause') -> List['ConstrainedAtom']:
    """For a given clause, return a list of all the constrained atoms in the clause"""
    constrained_atoms = []
    for literal in clause.unconstrained_clause.literals:
        constrained_atom = build_constrained_atom_from_literal(literal, clause)
        constrained_atoms.append(constrained_atom)
    return constrained_atoms


def build_constrained_atom_from_literal(literal, clause):
    """Build a constrained atom given a specific literal and its parent clause"""
    atom = literal.atom 
    positive_literal = Literal(atom, True)
    # constrained atoms subclass ConstrainedClause, so need to be built of an
    # unconstrained clause
    unconstrained_atom = UnconstrainedClause([positive_literal])
    constrained_atom = ConstrainedAtom(unconstrained_atom, clause.bound_vars, clause.cs)
    return constrained_atom


def have_same_predicate(c_atom1: 'ConstrainedAtom', c_atom2: 'ConstrainedAtom') -> bool:
    """For two constrained atoms, check if they have the same underlying predicate.

    NOTE: As currently written, the two predicates must be literally the same object, not just have the same
    name and arity"""
    predicate1 = c_atom1.unconstrained_clause.literals[0].atom.predicate
    predicate2 = c_atom2.unconstrained_clause.literals[0].atom.predicate
    return predicate1 == predicate2


def get_solutions(cs: 'ConstraintSet', variables: List['LogicalVariable']) -> List[List['Constant']]:
    """Get a list of solutions to the constraint set cs for specific variables.

    NOTE: This assumes that cs constains at least one InclusionConstraint.
    Also, we assume that there are no domain variables anywhere.
    Additionally, we only consider substituting in constants, not free variables."""
    variable_domains = get_variable_domains(cs, variables)

    # extract the equality and inequality constraints from the constraint set
    logical_constraints = get_all_logical_constraints(cs)
    sols: List[List['Constant']] = initiate_variable_recursion(variables, variable_domains, logical_constraints)
    return sols


def initiate_variable_recursion(variables: List['LogicalVariable'],
                                domains: List[Set['Constant']],
                                constraints: List['LogicalConstraint']) -> List[List['Constant']]:
    """Start recursion on each variable in turn to construct all the valid solutions to a constraint set"""
    # we pass through a dictionary that accumulates information about the variables
    variable_dict = {variable: {'domain': domain, 'equal_constants': set()} for variable, domain in zip(variables, domains)}
    return recurse(variables, variable_dict, constraints)


def recurse(remaining_variables: List['LogicalVariable'],
            variable_dict: Dict,
            constraints: List['LogicalConstraint']) -> List[List['Constant']]:
    """Recursively construct the solutions for each substitution to 'variable'"""

    # base case of recursion: no variables left
    if len(remaining_variables) == 0:
        return [[]]

    variable = remaining_variables[0]
    partial_sols: List[Tuple[List['Constant'], Dict]] = []
    for substitution in variable_dict[variable]['domain']: # iterate over allowed constants for current variable
        print(substitution)
        valid_substitution = True # to start with, assume the substitution will work out
        # create a copy of the dictionary to change for the new solution
        # NOTE: we don't copy the dict keys because these are LogicalVariables and we need them to be the same
        next_variable_dict = {key: deepcopy(value) for key, value in variable_dict.items()}
        constraints_involving_variable = get_relevant_logical_constraints(constraints, variable)

        # now follows a lot of messy logic to update the allowed constants for each variable
        for constraint in constraints_involving_variable:
            # first, deal with equality constraints that involve the current variable
            if isinstance(constraint, EqualityConstraint):
                if constraint.left_term == variable:
                    if isinstance(constraint.right_term, LogicalVariable) and constraint.right_term in remaining_variables: 
                        next_variable_dict[constraint.right_term]['equal_constants'].add(substitution)
                    elif isinstance(constraint.right_term, Constant) and constraint.right_term != substitution:
                        valid_substitution = False
                        break # exit this substitution and go back to substitution loop

                elif constraint.right_term == variable:
                    if isinstance(constraint.left_term, LogicalVariable) and constraint.left_term in remaining_variables: 
                        next_variable_dict[constraint.left_term]['equal_constants'].add(substitution)
                    elif isinstance(constraint.left_term, Constant) and constraint.left_term != substitution:
                        valid_substitution = False
                        break # exit this substitution and go back to substitution loop

            elif isinstance(constraint, InequalityConstraint):
                if constraint.left_term == variable:
                    if isinstance(constraint.right_term, LogicalVariable) and constraint.right_term in remaining_variables: 
                        next_variable_dict[constraint.right_term]['domain'].discard(substitution)
                    elif isinstance(constraint.right_term, Constant) and constraint.right_term == substitution:
                        valid_substitution = False
                        break # exit this substitution and go back to substitution loop

                elif constraint.right_term == variable:
                    if isinstance(constraint.left_term, LogicalVariable) and constraint.left_term in remaining_variables: 
                        next_variable_dict[constraint.left_term]['domain'].discard(substitution)
                    elif isinstance(constraint.left_term, Constant) and constraint.left_term == substitution:
                        valid_substitution = False
                        break # exit this substitution and go back to substitution loop
        # after all constraints are done, we check if this solution is still feasible
        if valid_substitution:
            partial_sols.append((substitution, next_variable_dict))

    # finally, we recurse on each partial solution
    all_solutions = []
    for sub, next_variable_dict in partial_sols:
        recursive_solutions = recurse(remaining_variables[1:], next_variable_dict, constraints)
        solutions = [[sub] + sol for sol in recursive_solutions]
        all_solutions += solutions
    # add all solutions together and return
    return all_solutions


def get_relevant_logical_constraints(constraints: List['LogicalConstraint'],
                                     variable: 'LogicalVariable') -> List['LogicalConstraint']:
    """From a list of logical constraints, get the ones that involve 'variable'"""
    return [constraint for constraint in constraints if (constraint.left_term == variable or constraint.right_term == variable)]


def get_variable_domains(cs: 'ConstraintSet', variables: List['LogicalVariable']) -> List[Set['Constant']]:
    """Find the allowed set of constants for each variable based on the set constraints
    
    NOTE: assumes there are no domain variables anywhere."""
    variable_domains = []
    for variable in variables:
        allowed_constants = get_variable_domain(cs, variable)
        variable_domains.append(allowed_constants)
    return variable_domains


def get_variable_domain(cs: 'ConstraintSet', variable: 'LogicalVariable') -> Set['Constant']:
    """Find the allowed set of constants for A SINGLE variable based on the set constraints
    
    NOTE: assumes there are no domain variables anywhere."""
    inclusion_constraints, notinclusion_constraints = get_relevant_set_constraints(cs.constraints, variable)
    assert(len(inclusion_constraints) > 0) # this should be true by assumption

    inclusion_constants = [ic.domain_term.constants for ic in inclusion_constraints]
    emptyset: Set['Constant'] = set() # a hack for type checking
    included_constants = emptyset.union(*inclusion_constants)

    # notinclusion constraints are optional
    if len(notinclusion_constraints) > 0: 
        notinclusion_constants = [nc.domain_term.constants for nc in notinclusion_constraints]
        notincluded_constants = emptyset.union(*notinclusion_constants)
    else:
        notincluded_constants = set()

    allowed_constants = included_constants - notincluded_constants
    return allowed_constants


def get_relevant_set_constraints(cs: List['Constraint'], variable: 'LogicalVariable') -> Tuple[List['InclusionConstraint'], List['NotInclusionConstraint']]:
    """Get lists of the inclusion and negated inclusion constraints from a constraint set for a specific variable"""
    inclusion_condition = lambda c: isinstance(c, InclusionConstraint) and c.logical_term == variable
    notinclusion_condition = lambda c: isinstance(c, NotInclusionConstraint) and c.logical_term == variable
    inclusion_constraints: List['InclusionConstraint'] = []
    notinclusion_constraints: List['NotInclusionConstraint']  = []
    for constraint in cs:
        if inclusion_condition(constraint):
            assert(isinstance(constraint, InclusionConstraint)) # assertion fixes type checking
            inclusion_constraints.append(constraint)
        elif notinclusion_condition(constraint):
            assert(isinstance(constraint, NotInclusionConstraint)) # assertion fixes type checking
            notinclusion_constraints.append(constraint)
    return inclusion_constraints, notinclusion_constraints


def get_all_logical_constraints(cs: 'ConstraintSet') -> List['LogicalConstraint']:
    """Get lists of the equality and inequality constraints from a constraint set"""
    return [constraint for constraint in cs.constraints if isinstance(constraint, LogicalConstraint)]


if __name__ == '__main__':
    preds = [Predicate('smokes', 2), Predicate('friends', 2), Predicate('fun', 1)]
    constants = [Constant('a'), Constant('b'), Constant('c')]
    variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]

    atoms = [Atom(preds[0], [variables[0], variables[1]]),
             Atom(preds[1], [variables[1], constants[0]]),
             Atom(preds[2], [constants[0]]),
             Atom(preds[1], [variables[1], variables[0]])
            ]

    literals = [Literal(atoms[0], True), Literal(atoms[1], False), Literal(atoms[2], False), Literal(atoms[3], True) ]
    domains = [SetOfConstants(constants), SetOfConstants(constants[:2])]

    constraint_set = ConstraintSet([EqualityConstraint(variables[0], constants[0]),
                      InequalityConstraint(variables[1], constants[1]),
                      InclusionConstraint(variables[0], domains[0]),
                      InclusionConstraint(variables[1], domains[0]),
                      NotInclusionConstraint(variables[1], domains[1])
                      ])
    print(constraint_set)
    print("get sol",get_solutions(constraint_set, variables[:2]))

    u_clause = UnconstrainedClause(literals)
    c_clause = ConstrainedClause(u_clause, variables[:2], constraint_set)

    c_atoms = get_constrained_atoms(c_clause)
    for c_atom in c_atoms:
        print(c_atom)

    print(have_same_predicate(c_atoms[1], c_atoms[2]))

