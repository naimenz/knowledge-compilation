"""File for messing around with closing substitution ideas"""

from kc.data_structures import *
from kc.util import *

from typing import Set


# a motivating example
X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
a, b, c = Constant('a'), Constant('b'), Constant('c')
D = SetOfConstants([a, b])
Universe = SetOfConstants([a, b, c])
p, q = Predicate('p', 1), Predicate('q', 1)
pX = Literal(Atom(p, [X]), True)
pY = Literal(Atom(p, [Y]), True)
qX = Literal(Atom(q, [X]), True)

XinD = InclusionConstraint(X, D)
YeqZ = EqualityConstraint(Y, Z)
XeqY = EqualityConstraint(X, Y)
Xeqa = EqualityConstraint(X, a)
Yeqa = EqualityConstraint(Y, a)
cs_gamma = ConstraintSet([XinD, ~YeqZ])
cs_gamma2 = ConstraintSet([XinD, YeqZ])
cs_A = ConstraintSet([XinD])

# sub = Substitution([(Z, Y)])
# new_cs_gamma = cs_gamma.apply_substitution(sub)
# print("hi",new_cs_gamma)

gamma = ConstrainedClause(UnconstrainedClause([pX, qX]), [X], cs_gamma)
A = ConstrainedAtom(UnconstrainedClause([pX]), [X], cs_A)

# print(gamma)
# print(get_free_logical_variables_in_clause(gamma))
# print(get_closing_substitutions(gamma, Universe))
cs1 = ConstraintSet([XinD, ~Xeqa])
cs2 = ConstraintSet([XinD, XeqY, Yeqa])

c_atom1 = ConstrainedAtom(UnconstrainedClause([pX]), [X], cs1)
c_atom2 = ConstrainedAtom(UnconstrainedClause([pX]), [X], cs2)
print(c_atom1,'\n',c_atom2)
indep = constrained_atoms_independent(c_atom1, c_atom2, Universe)
print(indep)
