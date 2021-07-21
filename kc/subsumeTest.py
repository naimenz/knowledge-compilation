# def subsumes(self, other){
# return split.is_applicable(self, other) && !split.is_applicable(other, self)
# }

def is_applicable(atom1: 'ConstrainedAtom', atom2: 'ConstrainedAtom'):
    mgu_eq_classes = atom1.get_constrained_atom_mgu_eq_classes(atom2)
    independent = True if mgu_eq_classes is None else False
    return not independent and atom2.does_not_subsume(atom1, mgu_eq_classes)

def does_not_subsume(self, other, mgu_eq_classes):
    mgu_substitution = mgu_eq_classes.to_substitution()
    this_atom = self.substitute(mgu_substitution)
    other_atom = other.substitute(mgu_substitution)
    if any(len(eq_class.constants) > 0
           and len(eq_class.variables.intersection(this_atom.bound_vars)) > 0
           for eq_class in mgu_eq_classes):
        return True
    elif any(len(eq_class.variables.intersection(this_atom.bound_vars)) >= 2
             for eq_class in mgu_eq_classes):
        return True
    elif any(inequality not in this_atom.get_constant_inequalities()
             for inequality in other_atom.get_constant_inequalities()):
        return True
    elif any(inequality not in this_atom.get_bound_variable_inequalities()
             and inequality.is_not_trivial_in(other_atom)
             for inequality in other_atom.get_bound_variable_inequalities()):
        return True
    else:
        return False


# TODO: put this in ConstrainedClause
# TODO: THINK ABOUT FREE VARIABLE CASE
def get_constant_inequalities(self) -> Set['NotInclusionConstraint']:
    ineq_constraints = set()
    for constraint in self.cs.set_constraints:
        if isinstance(constraint, NotInclusionConstraint):
            domain_term = constraint.domain_term
            if isinstance(domain_term, SetOfConstants) and domain_term.size() == 1:
                ineq_constraints.add(constraint)
    return ineq_constraints

# TODO: put this in ConstrainedClause
# TODO: THINK ABOUT FREE VARIABLE CASE
def get_bound_variable_inequalities(self) -> Set['InequalityConstraint']:
    ineq_constraints = set()
    for constraint in self.logical_constraints:
        if isinstance(constraint, InequalityConstraint):
            if constraint.left_term in self.bound_vars and constraint.right_term in self.bound_vars:
                ineq_constraints.add(constraint)
    return ineq_constraints

# TODO: put this in InequalityConstraint
def is_not_trivial(self, c_atom: 'ConstrainedAtom') -> bool:
    left_term_class = EquivalenceClass([self.left_term])
    right_term_class = EquivalenceClass([self.right_term])
    left_domain = left_term_class.get_shared_domain_from_cs(c_atom.cs)
    right_domain = right_term_class.get_shared_domain_from_cs(c_atom.cs)
    domain_terms_intersect = DomainTerm.intersection(left_domain, right_domain).size() > 0
    return domain_terms_intersect




















