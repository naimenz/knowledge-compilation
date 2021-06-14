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


def constrained_atoms_independent(c_atom1: 'ConstrainedAtom',
                                  c_atom2: 'ConstrainedAtom',
                                  universe: 'SetOfConstants'
                                  ) -> bool:
    """Are the constrained atoms c_atom1 and c_atom2 independent?
    NOTE: this can handle both formulas and sentences, so needs a universe"""
    c_atom1_sentences = build_sentences_from_open_formula(c_atom1, universe)
    c_atom2_sentences = build_sentences_from_open_formula(c_atom2, universe)
    # now check that all closing substitutions for the first are independent of all for the secone
    for sentence1 in c_atom1_sentences:
        for sentence2 in c_atom2_sentences:
            if not constrained_atoms_independent_sentences(sentence1, sentence2):
                return False
    return True

def build_sentences_from_open_formula(formula: C, universe: 'SetOfConstants') -> List[C]:
    """Build all possible sentences from a formula (either a sentence of non-sentence) and a given universe"""
    free_variables = get_free_logical_variables_in_clause(formula)
    # if not sentences, construct all the possible sentences from the formula
    if len(free_variables) > 0:
        closing_substitutions = get_closing_substitutions(formula, universe)
        sentences = [formula.apply_substitution(sub) for sub in closing_substitutions]
    else:
        sentences = [formula]
    return sentences


def constrained_atoms_independent_sentences(c_atom1: 'ConstrainedAtom', c_atom2: 'ConstrainedAtom') -> bool:
    """Are the SENTENCE constrained atoms c_atom1 and c_atom2 independent?
    This is checked by seeing if they have the same predicate and then whether they have the same groundings.
    Returns a boolean indicating whether they are.

    NOTE: currently assumes there are no free variables or domain variables"""
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


def is_satisfiable(constraint_set: 'ConstraintSet') -> bool:
    """Check if a constraint set is satisfiable by constructing solutions to it
    and seeing if there are any"""
    # first we check if any constraints contain a contradiction
    if constraint_set.contains_contradiction():
        return False
    # extract just the variables from each constraint
    logical_variables = get_logical_variables_from_cs(constraint_set)
    # now we find solutions for those variables and see if there are any
    solutions = get_solutions(constraint_set, list(logical_variables))
    return len(solutions) > 0


def get_logical_variables_from_cs(constraint_set: 'ConstraintSet') -> Set['LogicalVariable']:
    """Extract just the variables from each constraint in the constraint set"""
    logical_variables: Set['LogicalVariable'] = set()
    for constraint in constraint_set:
        for term in constraint.terms:
            if isinstance(term, LogicalVariable):
                logical_variables.add(term)
    return logical_variables


def get_free_logical_variables_in_clause(clause: 'ConstrainedClause') -> Set['LogicalVariable']:
    """Get only the free LOGICAL variables that appear in a clause.
    Free variables can appear in the constraint set OR the unconstrained clause."""

    all_variables = get_all_logical_variables_in_clause(clause)
    return all_variables - clause.bound_vars

def get_all_logical_variables_in_clause(clause: 'ConstrainedClause') -> Set['LogicalVariable']:
    """Get only the free LOGICAL variables that appear in a clause.
    Free variables can appear in the constraint set OR the unconstrained clause."""

    all_variables: Set['LogicalVariable'] = set()
    for literal in clause.unconstrained_clause.literals:
        all_variables.update(term for term in literal.atom.terms if isinstance(term, LogicalVariable))
    cs_logical_variables = get_logical_variables_from_cs(clause.cs)
    all_variables.update(cs_logical_variables)
    # return only the variables that are not quantified over
    return all_variables 


def get_closing_substitutions(clause: 'ConstrainedClause', universe: 'SetOfConstants') -> Set['Substitution']:
    """Get all satisfiable closing substitutions for a clause with free logical variables given a specific universe of constants."""
    free_logical_variables = get_free_logical_variables_in_clause(clause)
    # For a closing substitution, all free variables are substituted with constants or bound variables 
    bound_vars: Set['LogicalTerm'] = cast(Set['LogicalTerm'], clause.bound_vars) # hack for type checking
    universe_constants: Set['LogicalTerm'] =  cast(Set['LogicalTerm'], universe.constants) # hack for type checking
    substitution_targets: Set['LogicalTerm'] = bound_vars.union(universe_constants)
    # construct ALL possible closing substitutions first, then weed out the ones that are unsatisfiable
    right_hand_sides = product(substitution_targets, repeat=len(free_logical_variables))
    all_substitutions = set(Substitution(zip(free_logical_variables, rhs)) for rhs in right_hand_sides)
    valid_substitutions = set(sub for sub in all_substitutions if is_satisfiable(clause.cs.apply_substitution(sub)))
    return valid_substitutions


# DEBUG add variables argument
def get_constrained_atom_grounding(c_atom: 'ConstrainedAtom',
        variables: Optional[List['LogicalVariable']]=None) -> List['GroundAtom']:
    """Return a list of all ground atoms in the grounding of a constrained atom"""
    sols = get_solutions_to_constrained_atom(c_atom, variables)
    ground_atoms = [GroundAtom.build_from_atom_substitution(c_atom.atom, sol) for sol in sols]
    return ground_atoms

def constrained_atoms_subsumed(subsumer: 'ConstrainedAtom',
                                  subsumed: 'ConstrainedAtom',
                                  universe: 'SetOfConstants'
                                  ) -> bool:
    """Does the constrained atom subsumer subsume  the subsumed?
    NOTE: this can handle both formulas and sentences, so needs a universe"""
    print(f"original formulas: {subsumer},\n {subsumed}")
    # DEBUG
    subsumer_vars = get_all_logical_variables_in_clause(subsumer)
    subsumed_vars = get_all_logical_variables_in_clause(subsumed)

    subsumer_sentences = build_sentences_from_open_formula(subsumer, universe)
    subsumed_sentences = build_sentences_from_open_formula(subsumed, universe)
    # now check that all closing substitutions for the first are independent of all for the secone
    for subsumer_sentence in subsumer_sentences:
        for subsumed_sentence in subsumed_sentences:
            print(f"subsumer = {subsumer_sentence}\nsubsumed = {subsumed_sentence}")
            if not constrained_atoms_subsumed_sentences(subsumer_sentence, subsumed_sentence, list(subsumed_vars)):
                print("does NOT subsume\n")
                return False
            print("subsumes\n")
    return True


# DEBUG
def constrained_atoms_subsumed_sentences(subsumer: 'ConstrainedAtom', subsumed: 'ConstrainedAtom',
        variables: Optional[List['LogicalVariable']]=None) -> bool:
    """Check whether one SENTENCE constrained atom implies another.
    I.e. check whether the subsumer subsumes the subsumed"""
    # first, check if they have the same predicate, otherwise they cannot subsume
    if not have_same_predicate(subsumer, subsumed):
        return False
    
    subsumer_ground_atoms = get_constrained_atom_grounding(subsumer, variables)
    subsumed_ground_atoms = get_constrained_atom_grounding(subsumed, variables)

    # if there is a ground atom in the subsumed that does not appear in the subsumer, then
    # it is not subsumed
    for ground_atom in subsumed_ground_atoms:
        if not ground_atom in subsumer_ground_atoms:
            return False
    return True


def constrained_literals_subsumed(subsumer: 'UnitClause', subsumed: 'UnitClause', universe: 'SetOfConstants') -> bool:
    """Check whether one constrained literal implies another.
    This builds on constrained_atoms_subsumed with a check for polarity"""
    subsumer_constrained_atom = get_constrained_atoms(subsumer)[0] # only one literal
    subsumed_constrained_atom = get_constrained_atoms(subsumed)[0] # only one literal
    subsumed_as_atoms = constrained_atoms_subsumed(subsumer_constrained_atom, subsumed_constrained_atom, universe)
    polarity_same = (subsumer.literal.polarity == subsumed.literal.polarity)
    return subsumed_as_atoms and polarity_same


def constrained_clauses_subsumed(subsumer: 'ConstrainedClause', subsumed: 'ConstrainedClause', universe: 'SetOfConstants') -> bool:
    """Check whether one constrained clause implies another.
    This works by noticing that if any constrained literal of the subsumer implies any
    constrained literal of the subsumed, then it is subsumed."""
    subsumer_constrained_literals = get_constrained_literals(subsumer)
    subsumed_constrained_literals = get_constrained_literals(subsumed)

    for subsumer_literal in subsumer_constrained_literals:
        for subsumed_literal in subsumed_constrained_literals:
            if constrained_literals_subsumed(subsumer_literal, subsumed_literal, universe):
                return True
    return False


def constrained_atoms_not_independent_and_not_subsumed(subsumer: 'ConstrainedAtom', subsumed: 'ConstrainedAtom', universe: 'SetOfConstants') -> bool:
    return not constrained_atoms_independent(subsumer, subsumed, universe) and not constrained_atoms_subsumed(subsumer, subsumed, universe)


def constrained_clauses_independent_sentences(c_clause1: 'ConstrainedClause', c_clause2: 'ConstrainedClause') -> bool:
    """Are the constrained clauses c_clause1 and c_clause2 independent?
    This is checked by checking if every c-atom in the first clause is independent of every c-atom in the second.
    Returns a boolean indicating whether they are.

    NOTE: There may be ways to make this more efficient
    TODO: Replace this with non-sentence version"""
    c_atoms1 = get_constrained_atoms(c_clause1)
    c_atoms2 = get_constrained_atoms(c_clause2)
    for c_atom1 in c_atoms1:
        for c_atom2 in c_atoms2:
            if not constrained_atoms_independent_sentences(c_atom1, c_atom2):
                return False
    return True


# DEBUG added variables arg
def get_solutions_to_constrained_atom(c_atom: 'ConstrainedAtom',
                                      variables: Optional[List['LogicalVariable']]=None) -> Set['Substitution']:
    """Get the solutions to a constrained atom.
    NOTE: Assumes that the constrained atom contains no free or domain variables"""
    terms = c_atom.atom.terms

    if variables is None:
        # NOTE: we can only get solutions to variables, not constants
        variables = [term for term in terms if isinstance(term, LogicalVariable)]

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
                 ) -> Set['Substitution']:
    """Get a list of solutions to the constraint set cs for specific variables.
    Returns a list of lists of tuples of (variable, substitution) pairs

    NOTE: This assumes that cs constains at least one InclusionConstraint.
    Also, we assume that there are no domain variables anywhere.
    Additionally, we only consider substituting in constants, not free variables."""
    variable_domains = get_variable_domains(cs, variables)

    # extract the equality and inequality constraints from the constraint set
    logical_constraints = get_all_logical_constraints(cs)
    return set(initiate_variable_recursion(variables, variable_domains, logical_constraints))


def initiate_variable_recursion(variables: List['LogicalVariable'],
                                domains: List[Set['Constant']],
                                constraints: List['LogicalConstraint']
                                ) -> List['Substitution']:
    """Start recursion on each variable in turn to construct all the valid solutions to a constraint set"""
    # we pass through a dictionary that accumulates information about the variables
    variable_dict = {variable: {'domain': domain, 'substitution': None} for variable, domain in zip(variables, domains)}
    sols: List[Tuple[Dict, bool]] = recursively_construct_partial_sols(variables, variable_dict, constraints)
    return [build_substitution_from_variable_dict(sol) for sol, flag in sols if flag]


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
    
    constraints_involving_variable = get_relevant_logical_constraints(constraints, remaining_variables[0])
    if len(variable_dict[variable]['domain']) == 0:
        return [(dict(), False)] # no valid solutions for this variable, so the constraint set is unsatisfiable
    for potential_constant in variable_dict[variable]['domain']: # iterate over allowed constants for current variable
        next_variable_dict, valid_substitution = update_variable_dict_with_potential_constant(remaining_variables, variable_dict, constraints_involving_variable, potential_constant)
        if valid_substitution:
            partial_sols.append(next_variable_dict)

    # finally, we recurse on each partial solution
    return build_solutions(remaining_variables, constraints, partial_sols)


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
    # if there are still constraints, update those
    if len(constraints) > 0:
        next_variable_dict, valid_substitution = process_all_constraints(remaining_variables, variable_dict, constraints, potential_constant)
    # if not, then the substitution is valid and the variable dict stays the same
    else:
        valid_substitution = True
        next_variable_dict = copy_variable_dict(variable_dict)
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
                variable_dict[constraint.right_term]['domain'] = variable_dict[constraint.right_term]['domain'].intersection(set([potential_constant]))
            elif isinstance(constraint.right_term, Constant) and constraint.right_term != potential_constant:
                valid_substitution = False

        elif constraint.right_term == variable:
            if isinstance(constraint.left_term, LogicalVariable) and constraint.left_term in remaining_variables: 
                variable_dict[constraint.left_term]['domain'] = variable_dict[constraint.left_term]['domain'].intersection(set([potential_constant]))
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
    next_variable_dict = variable_dict
    for constraint in constraints:
        next_variable_dict, valid_substitution = process_constraint(remaining_variables, next_variable_dict, constraint, potential_constant)
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

    inclusion_constants = [set(ic.domain_term.constants) for ic in inclusion_constraints]
    included_constants = set.intersection(*inclusion_constants)

    # notinclusion constraints are optional
    if len(notinclusion_constraints) > 0: 
        notinclusion_constants = [set(nc.domain_term.constants) for nc in notinclusion_constraints]
        notincluded_constants = set.union(*notinclusion_constants)
    else:
        notincluded_constants = set()

    allowed_constants = included_constants - notincluded_constants
    return allowed_constants


def get_relevant_set_constraints(cs: FrozenSet['Constraint'], variable: 'LogicalVariable') -> Tuple[List['InclusionConstraint'], List['NotInclusionConstraint']]:
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


def build_substitution_from_variable_dict(variable_dict: Dict) -> 'Substitution':
    """The variable dict is the state that is modified during recursion.
    This function creates a Substitution object from that variable dict
    """
    variable_constant_pairs: List[VarTermPair] = [(var, subdict['substitution']) for var, subdict in variable_dict.items()]
    substitution = Substitution(variable_constant_pairs)
    return substitution
