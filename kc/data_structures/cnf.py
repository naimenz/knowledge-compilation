"""
Class for FO-CNF formulas.
"""

from kc.data_structures import EquivalenceClasses, UnconstrainedClause, ConstrainedClause, LogicalVariable, DomainVariable

from typing import List, Any, Iterable, Set
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import EquivalenceClass, Clause, Substitution, DomainTerm, Constant, SetOfConstants

class CNF:
    """
    A FOL-DC CNF.
    This consists of a set of (CONSTRAINED OR UNCONSTRAINED) clauses, which form a conjunction.
    """

    def __init__(self, clauses: Iterable['Clause'], shattered: bool = False) -> None:
        """Initialise with a set of clauses and (optionally) the shattering status of the cnf"""
        self.clauses = frozenset(clauses)

        u_clauses, c_clauses = set(), set()
        for clause in clauses:
            if isinstance(clause, UnconstrainedClause):
                u_clauses.add(clause)
            elif isinstance(clause, ConstrainedClause):
                c_clauses.add(clause)
        self.u_clauses = frozenset(u_clauses)
        self.c_clauses = frozenset(c_clauses)

        self.shattered = shattered # keep track of whether this cnf has undergone shattering

    def join(self, other: 'CNF') -> 'CNF':
        """Combine two CNFs into one."""
        return CNF(self.clauses.union(other.clauses))

    def substitute(self, substitution: 'Substitution') -> 'CNF':
        """Return a new CNF, the result of applying substitution to this CNF"""
        new_clauses = set(clause.substitute(substitution) for clause in self.clauses)
        return CNF(new_clauses)

    def get_unifying_classes(self) -> 'EquivalenceClasses':
        """Construct all unifying classes from this CNF
        and return them as EquivalenceClasse
        TODO: decide whether this should only consider bound variables (as per the definitions)
        or include free variables too (as per the examples and Forclift)"""
        initial_eq_class_pairs: List['EquivalenceClass'] = []
        for clause in self.clauses:
            for other_clause in self.clauses:
                for c_atom in clause.get_constrained_atoms():
                    for other_c_atom in other_clause.get_constrained_atoms():
                        eq_classes = c_atom.get_constrained_atom_mgu_eq_classes(other_c_atom)
                        if not eq_classes is None:
                            initial_eq_class_pairs += eq_classes
        initial_eq_classes = EquivalenceClasses(initial_eq_class_pairs)
        final_eq_classes = initial_eq_classes.make_self_consistent()
        return final_eq_classes

    def eq_class_has_one_variable(self, eq_class: 'EquivalenceClass') -> bool:
        """Determine whether a given root equivalence class has a single bound variable
        per clause or not.
        NOTE: We assume that this equivalence class is root in the cnf"""
        # first check if there are any unconstrained clauses, in which case we can't do this
        if len(self.u_clauses) > 0:
            return False

        for clause in self.c_clauses: 
            if len(eq_class.members.intersection(clause.bound_vars)) != 1:
                return False
        return True

    def eq_class_has_two_variables(self, eq_class: 'EquivalenceClass') -> bool:
        """Determine whether a given root equivalence class has two bound variables
        per clause or not.
        NOTE: We assume that this equivalence class is root in the cnf"""
        # first check if there are any unconstrained clauses, in which case we can't do this
        if len(self.u_clauses) > 0:
            return False

        for clause in self.c_clauses: 
            if len(eq_class.members.intersection(clause.bound_vars)) != 2:
                return False
        return True

    def get_logical_variables(self) -> Set['LogicalVariable']:
        """Return a set of all the logical variables that appear in the cnf"""
        logical_variables: Set['LogicalVariable'] = set()
        clause_logical_variables = [clause.all_variables for clause in self.clauses]
        logical_variables.update(*clause_logical_variables)
        return logical_variables

    def get_free_logical_variables(self) -> Set['LogicalVariable']:
        """Return a set of all FREE logical variables in the cnf
        NOTE: kind of assumes each variable appears in one clause"""
        all_variables = self.get_logical_variables()
        bound_variables: Set['LogicalVariable']  = set() # hack for type checking
        clause_bound_variables = [cc.bound_vars for cc in self.c_clauses]
        bound_variables.update(*clause_bound_variables)
        return all_variables - bound_variables

    def get_constants(self) -> Set['Constant']:
        """Get all constants that appear as logical terms (i.e. not in domain terms) in the cnf"""
        constants: Set['Constant'] = set()
        clause_constants = [clause.constants for clause in self.clauses]
        constants.update(*clause_constants)
        return constants

    def get_domain_terms(self) -> Set['DomainTerm']:
        """Get all the domain terms that appear in the cnf"""
        domain_terms = set()
        for clause in self.c_clauses: # u_clauses have no domain terms
            for constraint in clause.cs.set_constraints:
                domain_terms.add(constraint.domain_term)
        return domain_terms

    def get_domain_variables(self) -> Set['DomainVariable']:
        """Get all the domain VARIABLES that appear in the cnf"""
        domain_variables = set()
        for clause in self.c_clauses: # u_clauses have no domain terms
            for constraint in clause.cs.set_constraints:
                if isinstance(constraint.domain_term, DomainVariable):
                    domain_variables.add(constraint.domain_term)
        return domain_variables
        
    def get_new_logical_variable(self, symbol: str) -> 'LogicalVariable':
        """Return a logical variable that does not appear in the theory.
        To make it unique, take the symbol and keep adding underscores"""
        new_variable_string = symbol
        new_variable = LogicalVariable(new_variable_string)
        logical_variables = self.get_logical_variables()
        while new_variable in logical_variables:
            new_variable_string += '_'
            new_variable = LogicalVariable(new_variable_string)
        return new_variable

    def get_new_domain_variable(self,
                                symbol: str,
                                parent_domain: 'ProperDomain',
                                excluded_constants: Set['Constant']
                                ) -> 'DomainVariable':
        """Return a logical variable that does not appear in the theory.
        To make it unique, take the symbol and keep adding underscores"""
        new_variable_string = symbol
        domain_variable_symbols = [d.symbol for d in self.get_domain_variables()]
        while new_variable_string in domain_variable_symbols:
            new_variable_string += '_'
        new_variable = DomainVariable(new_variable_string, parent_domain, excluded_constants)
        return new_variable

    def __eq__(self, other: Any) -> bool:
        """Two CNFs are equal if they have the same clauses"""
        return isinstance(other, CNF) and self.clauses == other.clauses

    def __hash__(self) -> int:
        return hash(self.clauses)

    def __str__(self) -> str:
        clause_strs = [f'({str(clause)})' for clause in self.clauses]
        return '\nAND\n'.join(clause_strs)

    def __repr__(self) -> str:
        return self.__str__()

