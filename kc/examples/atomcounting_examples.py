from kc.data_structures import *
from kc.compiler import *

X = LogicalVariable('X')
Y = LogicalVariable('Y')
Z = LogicalVariable('Z')
X1 = LogicalVariable('X1')
Y1 = LogicalVariable('Y1')
Z1 = LogicalVariable('Z1')
X2 = LogicalVariable('X2')


alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
fun = Predicate('fun', 1)

friendsXY = Literal(Atom(friends, [X, Y]))
funX = Literal(Atom(fun, [X]))

friendsY1X1 = Literal(Atom(friends, [Y1, X1]))
funX1 = Literal(Atom(fun, [X1]))

People = RootDomain([alice, bob, charlie], 'People')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)

cs1 = ConstraintSet([XinPeople, YinPeople])
cs2 = ConstraintSet([X1inPeople, Y1inPeople])

clause1 = ConstrainedClause([funX, ~friendsXY], [X, Y], cs1)
clause2 = ConstrainedClause([funX1, ~friendsY1X1], [X1, Y1], cs2)

D = DomainVariable('D', People)
XinD = InclusionConstraint(X, D)
X1inD = InclusionConstraint(X1, D)
u1 = ConstrainedAtom([funX], [X], ConstraintSet([XinPeople, XinD]))
u2 = UnitClause([~funX1], [X1], ConstraintSet([X1inPeople, ~X1inD]))
# print(f'{u1=}')
# for ca in clause1.get_constrained_atoms():
#     print(f"================= {ca} =============")
#     print(ca.needs_splitting(u1))
cnf = CNF([clause1, clause2])
cnf.shattered = True  # hack so I can test AC directly like in the PhD, when really should shatter first
compiler = Compiler()
nnf = compiler.compile(cnf)
