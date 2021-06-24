"""Examples from the data_structures files all in one place and not imported"""

from kc.data_structures import *
from typing import List, Any

# logicalterms.py
print('logicalterms.py')
c = Constant('a')
print(c.value)
print(c)

x = LogicalVariable('X')
print(x.symbol)
print(x)
print([x])


# literals.py
print('literals.py')
pred = Predicate('smokes', 4)
pred2 = Predicate('smokes', 4)
print(pred)
print(pred == pred2)

terms: List[Any] = [Constant('a'), LogicalVariable('X'), Constant('b'), LogicalVariable('Y')]
atom = Atom(pred, terms)
print(atom)

ground_atom = Atom(pred, [terms[0], terms[2], terms[0], terms[0]])
print("ground atom",ground_atom)

substitution = Substitution([(terms[1], terms[0]), (terms[3], terms[2])])
built_atom = GroundAtom.build_from_atom_substitution(atom, substitution)
print("built atom", built_atom)

literal = Literal(atom, False)
print(literal)


# domainterms.py
print('domainterms.py')
c1 = Constant('a')
c2 = Constant('b')
c3 = Constant('c')

v1 = LogicalVariable('X')

D = SetOfConstants([c1, c2, c3])
print(D)

dom_var1 = DomainVariable('D')
dom_var2 = DomainVariable('D')
dom_var3 = dom_var1
print(dom_var1 == dom_var2)
print(dom_var1 == dom_var3)


# constraints.py
print('constraints.py')
v1 = LogicalVariable('X')
c1 = Constant('bob')
c2 = Constant('a')

eq_constraint = EqualityConstraint(v1, c2)
print(eq_constraint)
eq_constraint1 = EqualityConstraint(c2, v1)
print(eq_constraint1)
ineq_constraint = InequalityConstraint(v1, c2)
print(ineq_constraint)
print(eq_constraint == ineq_constraint)
print(eq_constraint == eq_constraint1)

d = SetOfConstants([c1, c2])
in_constraint = InclusionConstraint(c1, d)
print(in_constraint)

notin_constraint = NotInclusionConstraint(c1, d)
print(notin_constraint)
print("negated",~notin_constraint)

cs = ConstraintSet([eq_constraint, ineq_constraint, in_constraint, notin_constraint])
print(cs)



# clauses.py
print('clauses.py')
pred1 = Predicate('smokes', 1)
pred2 = Predicate('friends', 2)

c1 = Constant('bob')
v1 = LogicalVariable('X')
atom1 = Atom(pred1, [c1])
atom2 = Atom(pred2, [c1, v1])

literal1 = Literal(atom1, True)
literal2 = Literal(atom2, False)

uclause = UnconstrainedClause([literal1, literal2])

# v1 = LogicalVariable('X')
# c1 = Constant('bob')
c2 = Constant('a')
d = SetOfConstants([c1, c2])

eq_constraint = EqualityConstraint(v1, c2)
ineq_constraint = InequalityConstraint(v1, c2)
in_constraint = InclusionConstraint(c1, d)
notin_constraint = NotInclusionConstraint(c1, d)

cs = ConstraintSet([eq_constraint, ineq_constraint, in_constraint, notin_constraint])

v2 = LogicalVariable('Y')
clause = ConstrainedClause(uclause, [v1, v2], cs)
print(clause)


# cnf.py
print('cnf.py')
pred1 = Predicate('smokes', 1)
pred2 = Predicate('friends', 2)

c1 = Constant('bob')
v1 = LogicalVariable('X')
atom1 = Atom(pred1, [c1])
atom2 = Atom(pred2, [c1, v1])

literal1 = Literal(atom1, True)
literal2 = Literal(atom2, False)

uclause = UnconstrainedClause([literal1, literal2])

# v1 = LogicalVariable('X')
# c1 = Constant('bob')
c2 = Constant('a')
d = SetOfConstants([c1, c2])

eq_constraint = EqualityConstraint(v1, c2)
ineq_constraint = InequalityConstraint(v1, c2)
in_constraint = InclusionConstraint(c1, d)
notin_constraint = NotInclusionConstraint(c1, d)

cs1 = ConstraintSet([eq_constraint, in_constraint])
cs2 = ConstraintSet([ineq_constraint, notin_constraint])

v2 = LogicalVariable('Y')
v3 = LogicalVariable('Z')
clause1 = ConstrainedClause(uclause, [v1, v2], cs1)
clause2 = ConstrainedClause(uclause, [v2, v3], cs2)

cnf = CNF([clause1, clause2])
cnf2 = CNF([clause1, clause1])
print(cnf)
print(cnf == cnf)
print(cnf == cnf2)


# substitutions.py
print('substitutions.py')
variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]
sub_constants: List['LogicalTerm'] = [Constant('a'), Constant('b'), Constant('a')]

pairs = [(v, c) for v, c in zip(variables, sub_constants)]
substitution = Substitution(pairs)
print(substitution)
print(substitution[variables[0]])
print(substitution == substitution)
sub2 = Substitution(pairs[1::-1])
print(sub2)
print(substitution == sub2)

for abab in substitution.mappings():
    print(abab)

for abab in substitution:
    print(abab)

# equivalenceclasses.py
print('equivalenceclasses.py')
a, b, c = Constant('a'), Constant('b'), Constant('c')
X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z') 
eq_class1 = EquivalenceClass([a, b])
eq_class2 = EquivalenceClass([X, Y])
print(eq_class1.is_inconsistent)
print(eq_class2.is_inconsistent)
print(eq_class1.join(eq_class2).is_consistent)
