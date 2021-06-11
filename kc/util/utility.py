"""
This file contains various useful functions for working with the formulas defined in data_structures
"""

from kc.data_structures.cnf import *
from kc.data_structures.literals import *
from kc.data_structures.constraints import *
from kc.data_structures.logicalterms import *
from kc.data_structures.substitutions import *

from copy import deepcopy
from itertools import chain

from typing import List, Tuple, Set, Dict, cast

# defining type alias to simplify type hinting
VarTermPair = Tuple['LogicalVariable', 'LogicalTerm']


def get_constrained_atom_mgu(c_atom1: 'ConstrainedAtom', c_atom2: 'ConstrainedAtom') -> 'Substitution':
    """Get the mgu (most general unifier) between two constrained atoms, and return it as a substitution.

    NOTE: This assumes that the two c_atoms are not independent, which should have been checked earlier."""
    # first, try to unify the unconstrained atoms
    u_atom1 = c_atom1.atom
    u_atom2 = c_atom2.atom
    var_term_pair: VarTermPair
    var_term_pairs: List[VarTermPair] = []
    # TODO: Really think through how I should construct this substitution
    for term1, term2 in zip(u_atom1.terms, u_atom2.terms):
        if term1 != term2:
            # since we assumed the two aren't independent, at least one of these terms must be a variable
            if isinstance(term1, LogicalVariable):
                # check that we aren't already substituting these variables
                if any(term1 in pair for pair in var_term_pairs):
                    continue
                if isinstance(term2, LogicalVariable):
                    if any(term2 in pair for pair in var_term_pairs):
                        continue
                var_term_pair = (term1, term2)
            elif isinstance(term1, Constant) and isinstance(term2, LogicalVariable):
                if any(term1 in pair for pair in var_term_pairs):
                    var_term_pair = (term2, term1)
            else:
                raise ValueError('c_atoms {c_atom1}, {c_atom2} should not be independent')
            var_term_pairs.append(var_term_pair)


def get_unconstrained_atom_mgu(u_atom1: 'Atom', u_atom2: 'Atom') -> 'Substitution':
    """Get the mgu between two unconstrained atoms and return it as a substitution.

    NOTE: This assumes that the two u_atoms are not independent"""
    var_term_pair: VarTermPair
    var_term_pairs: List[VarTermPair] = []
    # TODO: Really think through how I should construct this substitution
    for term1, term2 in zip(u_atom1.terms, u_atom2.terms):
        if term1 != term2:
            # since we assumed the two aren't independent, at least one of these terms must be a variable
            if isinstance(term1, LogicalVariable):
                # check that we aren't already substituting the left variable
                if variable_already_substituted_incompatibly(term1, term2, var_term_pairs):
                    pass
                var_term_pair = (term1, term2)
            elif isinstance(term1, Constant) and isinstance(term2, LogicalVariable):
                var_term_pair = (term2, term1)
            else:
                raise ValueError('c_atoms {c_atom1}, {c_atom2} should not be independent')
            var_term_pairs.append(var_term_pair)

    print(var_term_pairs)
    print(Substitution(var_term_pairs))
    return simplify_substitution(Substitution(var_term_pairs))

def variable_already_substituted_incompatibly(variable: 'LogicalVariable',
                                 term: 'LogicalTerm',
                                 var_term_pairs: List[VarTermPair]
                                 ) -> bool:
    """Given a variable, a term and a list of tuples for substitutions, return True if the variable is 
    already being substituted to some term that is not 'term'"""
    for pair in var_term_pairs:
        if pair[0] == variable:
            return True
    return False


def simplify_substitution(substitution: 'Substitution') -> 'Substitution':
    """Sometimes substitutions contain redundancies that can be removed with transitivity.
    This function finds those redundancies and returns an equivalent substitution with them removed.

    NOTE: This assumes that substitutions are applied 'iteratively' until the formula doesn't change any more.
    """
    # the key here is to identify variables that appear on the left of some mapping and the right of some other mapping
    overlap = find_repeated_variables_in_substitution(substitution)
    while len(overlap) > 0:
        # grab a random variable that appears on both sides
        current_var = overlap.pop() 
        substitution = remove_transitivity_redundant_mappings(substitution, current_var)
        substitution = remove_equality_redundant_mappings(substitution)
        overlap = find_repeated_variables_in_substitution(substitution)
    # if there is no overlap, we are done and return straight away
    return substitution

def find_repeated_variables_in_substitution(substitution: 'Substitution') -> Set['LogicalVariable']:
    """Return a set of variables that appear on the left of some mapping and on the right of some other mapping
    in the substitution"""
    left_vars = set()
    right_vars = set()
    for mapping in substitution.mappings():
        left_var = mapping[0]
        right_term = mapping[1]
        left_vars.add(left_var)
        if isinstance(right_term, LogicalVariable):
            right_vars.add(right_term)
    overlap = left_vars.intersection(right_vars)
    return overlap

def remove_transitivity_redundant_mappings(substitution: 'Substitution', variable: 'LogicalVariable') -> 'Substitution':
    """Return a substitution with transitivity redundant (e.g. X -> Y, Y -> Z)  mappings involving a specific variable removed"""
    # find mappings where variable is on the right and one where it is on the left
    # and save the rest to a new list
    saved_mappings = []
    right_mappings = []
    for mapping in substitution.mappings():
        if mapping[0] == variable:
            found_transitive_redundancy = True
            left_mapping = mapping
            saved_mappings.append(left_mapping)
        elif mapping[1] == variable:
            right_mappings.append(mapping)
        else:
            saved_mappings.append(mapping)
    # add a new mapping for each right mapping that we remove,
    # simplifying by transitivity
    new_mappings = [(right_mapping[0], left_mapping[1]) for right_mapping in right_mappings]
    # build a new simplified substitution
    substitution = Substitution(saved_mappings + new_mappings)
    return substitution

def remove_equality_redundant_mappings(substitution: 'Substitution') -> 'Substitution':
    """Return a substitution with all equality redundant (e.g. X = X) mapping removed"""
    saved_mappings = []
    for mapping in substitution.mappings():
        if mapping[0] != mapping[1]:
            saved_mappings.append(mapping)
    return Substitution(saved_mappings)




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


def constrained_atoms_independent(c_atom1: 'ConstrainedAtom', c_atom2: 'ConstrainedAtom') -> bool:
    """Are the constrained atoms c_atom1 and c_atom2 independent?
    This is checked by seeing if they have the same predicate and then whether they have the same groundings.
    Returns a boolean indicating whether they are.

    NOTE: currently assumes there are no constants in the arguments to either atom"""
    # first, check if they have the same predicate, otherwise they are definitely independent
    if not have_same_predicate(c_atom1, c_atom2):
        return True
    
    ground_atoms1 = get_constrained_atom_grounding(c_atom1)
    ground_atoms2 = get_constrained_atom_grounding(c_atom2)

    # if they share a ground atom, then they are not independent
    for ground_atom in ground_atoms1:
        if ground_atom in ground_atoms2:
            return False
    return True

def get_constrained_atom_grounding(c_atom: 'ConstrainedAtom') -> List['GroundAtom']:
    """Return a list of all ground atoms in the grounding of a constrained atom"""
    sols = get_solutions_to_constrained_atom(c_atom)
    ground_atoms = [GroundAtom.build_from_atom_substitution(c_atom.atom, sol) for sol in sols]

    return ground_atoms

def constrained_atoms_subsumed(subsumer: 'ConstrainedAtom', subsumed: 'ConstrainedAtom') -> bool:
    """Check whether one constrained atom implies another.
    I.e. check whether the subsumer subsumes the subsumed"""
    # first, check if they have the same predicate, otherwise they cannot subsume
    if not have_same_predicate(subsumer, subsumed):
        return False
    
    subsumer_ground_atoms = get_constrained_atom_grounding(subsumer)
    subsumed_ground_atoms = get_constrained_atom_grounding(subsumed)

    # if there is a ground atom in the subsumed that does not appear in the subsumer, then
    # it is not subsumed
    for ground_atom in subsumed_ground_atoms:
        if not ground_atom in subsumer_ground_atoms:
            return False
    return True


def constrained_clauses_independent(c_clause1: 'ConstrainedClause', c_clause2: 'ConstrainedClause') -> bool:
    """Are the constrained clauses c_clause1 and c_clause2 independent?
    This is checked by checking if every c-atom in the first clause is independent of every c-atom in the second.
    Returns a boolean indicating whether they are.

    NOTE: There may be ways to make this more efficient"""
    c_atoms1 = get_constrained_atoms(c_clause1)
    c_atoms2 = get_constrained_atoms(c_clause2)
    for c_atom1 in c_atoms1:
        for c_atom2 in c_atoms2:
            if not constrained_atoms_independent(c_atom1, c_atom2):
                return False
    return True


def get_solutions_to_constrained_atom(c_atom: 'ConstrainedAtom') -> List['Substitution']:
    """NOTE: for now we assume that all the arguments to the constrained atom are variables"""
    terms = c_atom.unconstrained_clause.literals[0].atom.terms

    # TODO: remove this, it's just a hack for type-checking while I can't handle Constants
    assert(all(isinstance(term, LogicalVariable) for term in terms)) # just ensuring that we only pass variables
    variables: List['LogicalVariable'] = cast(List['LogicalVariable'], terms)

    # now we get the solutions for its constraint set with its variables
    return get_solutions(c_atom.cs, variables)


def have_same_predicate(c_atom1: 'ConstrainedAtom', c_atom2: 'ConstrainedAtom') -> bool:
    """For two constrained atoms, check if they have the same underlying predicate.

    NOTE: As currently written, the two predicates must be literally the same object, not just have the same
    name and arity"""
    predicate1 = c_atom1.atom.predicate
    predicate2 = c_atom2.atom.predicate
    return predicate1 == predicate2


def get_solutions(cs: 'ConstraintSet',
                  variables: List['LogicalVariable']
                 ) -> List['Substitution']:
    """Get a list of solutions to the constraint set cs for specific variables.
    Returns a list of lists of tuples of (variable, substitution) pairs

    NOTE: This assumes that cs constains at least one InclusionConstraint.
    Also, we assume that there are no domain variables anywhere.
    Additionally, we only consider substituting in constants, not free variables."""
    variable_domains = get_variable_domains(cs, variables)

    # extract the equality and inequality constraints from the constraint set
    logical_constraints = get_all_logical_constraints(cs)
    return initiate_variable_recursion(variables, variable_domains, logical_constraints)


def initiate_variable_recursion(variables: List['LogicalVariable'],
                                domains: List[Set['Constant']],
                                constraints: List['LogicalConstraint']
                                ) -> List['Substitution']:
    """Start recursion on each variable in turn to construct all the valid solutions to a constraint set"""
    # we pass through a dictionary that accumulates information about the variables
    variable_dict = {variable: {'domain': domain, 'equal_constants': set(), 'substitution': None} for variable, domain in zip(variables, domains)}
    sols: List[Tuple[Dict, bool]] = recursively_construct_partial_sols(variables, variable_dict, constraints)
    return [build_substitution_from_variable_dict(sol) for sol, flag in sols]


def build_substitution_from_variable_dict(variable_dict: Dict) -> 'Substitution':
    """The variable dict is the state that is modified during recursion.
    This function creates a Substitution object from that variable dict
    """
    variable_constant_pairs: List[VarTermPair] = [(var, subdict['substitution']) for var, subdict in variable_dict.items()]
    substitution = Substitution(variable_constant_pairs)
    return substitution


def recursively_construct_partial_sols(remaining_variables: List['LogicalVariable'],
            variable_dict: Dict,
            constraints: List['LogicalConstraint']
            ) -> List[Tuple[Dict, bool]]:
    """Recursively construct the solutions for each substitution to 'variable'
    Returns a list of tuples of a state dict and a flag for satisfiability or unsatisfiability"""

    # base case of recursion: no variables left
    if len(remaining_variables) == 0:
        return [(variable_dict, True)]

    variable = remaining_variables[0]
    partial_sols: List[Dict] = []
    
    # before considering substitutions, we need to process equality constraints on the variable
    variable_dict, satisfiable = process_equalities(variable, variable_dict)
    if not satisfiable:
        return [(dict(), False)]

    constraints_involving_variable = get_relevant_logical_constraints(constraints, remaining_variables[0])
    for potential_constant in variable_dict[variable]['domain']: # iterate over allowed constants for current variable
        next_variable_dict, valid_substitution = update_variable_dict_with_potential_constant(remaining_variables, variable_dict, constraints_involving_variable, potential_constant)
        if valid_substitution:
            partial_sols.append(next_variable_dict)

    # finally, we recurse on each partial solution
    return build_solutions(remaining_variables, constraints, partial_sols)


def process_equalities(variable: 'LogicalVariable', variable_dict: Dict) -> Tuple[Dict, bool]:
    """This function is called before calculating solutions for 'variable'.
    It propagates the repurcussions of equality constraints on the allowed values of 'variable'."""
    variable_dict = copy_variable_dict(variable_dict)

    if len(variable_dict[variable]['equal_constants']) > 1:
        return dict(), False
    if len(variable_dict[variable]['equal_constants']) == 1:
        variable_dict[variable]['domain'] = variable_dict[variable]['domain'].intersection(variable_dict[variable]['equal_constants'])
    return variable_dict, True


def copy_variable_dict(variable_dict: Dict) -> Dict:
    """Copies a variable dict to be mutated for the next recursion step.
    This is somewhere between a shallow and a deep copy -- the LogicalVariable keys stay the same
    but the sets of Constants are deep copied"""
    return {key: deepcopy(value) for key, value in variable_dict.items()}


def update_variable_dict_with_potential_constant(remaining_variables: List['LogicalVariable'],
                                  variable_dict: Dict,
                                  constraints: List['LogicalConstraint'],
                                  potential_constant: 'Constant') -> Tuple[Dict, bool]:
    """Update a variable dict with constraint information for a given potential constant.
    Returns an updated variable dict and a flag to say whether the constraint set is satisfiable with this potential constant"""
    next_variable_dict, valid_substitution = process_all_constraints(remaining_variables, variable_dict, constraints, potential_constant)
    if valid_substitution: # only update the dict if the substitution is valid (otherwise next_variable_dict isn't used)
        next_variable_dict[remaining_variables[0]]['substitution'] = potential_constant
    return next_variable_dict, valid_substitution


def process_constraint(remaining_variables: List['LogicalVariable'], variable_dict: Dict, constraint: 'LogicalConstraint', potential_constant: 'Constant') -> Tuple[Dict, bool]:
    """Update a variable dict with the implications of a specific constraint for a specific substitution
    and return the new variable dict and a flag indicating if the constraint is satisfiable"""
    variable = remaining_variables[0]
    valid_substitution = True

    variable_dict = copy_variable_dict(variable_dict)
    # TODO: refactor this horrible nested if-else structure
    if isinstance(constraint, EqualityConstraint):
        if constraint.left_term == variable:
            if isinstance(constraint.right_term, LogicalVariable) and constraint.right_term in remaining_variables: 
                variable_dict[constraint.right_term]['equal_constants'].add(potential_constant)
            elif isinstance(constraint.right_term, Constant) and constraint.right_term != potential_constant:
                valid_substitution = False

        elif constraint.right_term == variable:
            if isinstance(constraint.left_term, LogicalVariable) and constraint.left_term in remaining_variables: 
                variable_dict[constraint.left_term]['equal_constants'].add(potential_constant)
            elif isinstance(constraint.left_term, Constant) and constraint.left_term != potential_constant:
                valid_substitution = False

    elif isinstance(constraint, InequalityConstraint):
        if constraint.left_term == variable:
            if isinstance(constraint.right_term, LogicalVariable) and constraint.right_term in remaining_variables: 
                variable_dict[constraint.right_term]['domain'].discard(potential_constant)
            elif isinstance(constraint.right_term, Constant) and constraint.right_term == potential_constant:
                valid_substitution = False

        elif constraint.right_term == variable:
            if isinstance(constraint.left_term, LogicalVariable) and constraint.left_term in remaining_variables: 
                variable_dict[constraint.left_term]['domain'].discard(potential_constant)
            elif isinstance(constraint.left_term, Constant) and constraint.left_term == potential_constant:
                valid_substitution = False
    return variable_dict, valid_substitution


def process_all_constraints(remaining_variables: List['LogicalVariable'],
                            variable_dict: Dict,
                            constraints: List['LogicalConstraint'],
                            potential_constant: 'Constant') -> Tuple[Dict, bool]:
    """Process a batch of constraints and return an updated variable dict"""
    for constraint in constraints:
        next_variable_dict, valid_substitution = process_constraint(remaining_variables, variable_dict, constraint, potential_constant)
        if not valid_substitution: # save time by moving on to next substitution
            return dict(), False
    return next_variable_dict, True


def build_solutions(remaining_variables: List['LogicalVariable'],
                    constraints: List['LogicalConstraint'],
                    partial_sols: List[Dict]) -> List[Tuple[Dict, bool]]:
    """Construct solutions recursively using the valid substitutions found and the information propagated in the variable dicts"""
    all_solutions = []
    for next_variable_dict in partial_sols:
        recursive_solutions = recursively_construct_partial_sols(remaining_variables[1:], next_variable_dict, constraints)
        all_solutions += recursive_solutions
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
    # preds = [Predicate('smokes', 2), Predicate('friends', 2), Predicate('fun', 1)]
    # constants: List['Constant'] = [Constant('a'), Constant('b'), Constant('c')]
    # variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]

    # atoms = [Atom(preds[0], [variables[0], variables[1]]),
    #          Atom(preds[1], [variables[1], constants[0]]),
    #          Atom(preds[2], [constants[0]]),
    #          Atom(preds[1], [variables[0], variables[1]]),
    #          Atom(preds[1], [variables[1], variables[0]]),
    #         ]

    # literals = [Literal(atoms[0], True), Literal(atoms[1], False), Literal(atoms[2], False), Literal(atoms[3], True), Literal(atoms[4], True) ]
    # domains = [SetOfConstants(constants), SetOfConstants(constants[:2])]

    # constraints = [
    #                   # EqualityConstraint(variables[0], constants[0]),
    #                   InequalityConstraint(variables[1], constants[1]),
    #                   InequalityConstraint(variables[0], variables[1]),
    #                   InclusionConstraint(variables[0], domains[0]),
    #                   InclusionConstraint(variables[1], domains[0]),
    #                   # NotInclusionConstraint(variables[1], domains[1])
    #                   ]
    # constraint_set = ConstraintSet(constraints)
    # # print(constraint_set)
    # # print("get sol",get_solutions(constraint_set, variables[:2]))


    # constraints2 =[
    #                   # EqualityConstraint(variables[0], constants[0]),
    #                   InequalityConstraint(variables[1], constants[1]),
    #                   InequalityConstraint(variables[0], variables[1]),
    #                   InclusionConstraint(variables[0], domains[0]),
    #                   InclusionConstraint(variables[1], domains[0]),
    #                   # NotInclusionConstraint(variables[1], domains[1])
    #                   ]

    # constraint_set2 = ConstraintSet(constraints2)
    # # print(constraint_set2)
    # # print("get sol",get_solutions(constraint_set2, variables))


    # u_clause = UnconstrainedClause(literals[0:1] + literals[3:])
    # c_clause1 = ConstrainedClause(u_clause, variables[:2], constraint_set)
    # c_clause2 = ConstrainedClause(u_clause, variables[:2], constraint_set2)

    # c_atom1 = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], constraint_set)
    # c_atom2 = ConstrainedAtom(UnconstrainedClause([literals[-1]]), variables[:2], constraint_set2)
    # print("C_ATOMs:\n", c_atom1,'\n', c_atom2)
    # # print("ARGUMENTS:", get_solutions_to_constrained_atom(c_atom1))
    # print("Are the c-atoms independent?",constrained_atoms_independent(c_atom1, c_atom2))
    # print("Are the c-atoms subsumed?",constrained_atoms_subsumed(c_atom1, c_atom2))
    # print(get_constrained_atom_grounding(c_atom1))
    # print(get_constrained_atom_grounding(c_atom2))

    # # print(have_same_predicate(c_atoms[1], c_atoms[2]))
    # print("C_Clauses:\n", c_clause1,'\n',c_clause2)
    # print("Are the c-clauses independent?", constrained_clauses_independent(c_clause1, c_clause2))

    # subsumed = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], constraint_set)
    # subsumer = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], ConstraintSet(constraints[1:]))
    # print("Are the c-atoms subsumed?",constrained_atoms_subsumed(subsumer, subsumed))
    # print("Are the c-atoms subsumed?",constrained_atoms_subsumed(subsumed, subsumer))

    p_pred = Predicate('p', 3)
    a, b = Constant('a'), Constant('b')
    X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
    u_atom1 = Atom(p_pred, [X, X, a])
    u_atom2 = Atom(p_pred, [Y, Z, Z])
    print("MGU print")
    sub = get_unconstrained_atom_mgu(u_atom1, u_atom2)
    print(sub)

    # print(get_constrained_atom_mgu(subsumer, subsumed))
