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
dislikes = Predicate('dislikes', 2)
likes = Predicate('likes', 2)

friendsXY = Literal(Atom(friends, [X, Y]), True)
dislikesXY = Literal(Atom(dislikes, [X, Y]), True)

friendsX1Y1 = Literal(Atom(friends, [X1, Y1]), True)
likesX1Y1 = Literal(Atom(likes, [X1, Y1]), True)

friendsX2X2 = Literal(Atom(friends, [X2, X2]), True)

People = SetOfConstants([alice, bob, charlie])

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)

X2inPeople = InclusionConstraint(X2, People)

XeqY = EqualityConstraint(X, Y)
X2eqY = EqualityConstraint(X2, Y)
X2eqX = EqualityConstraint(X2, X)

cs = ConstraintSet([XinPeople, YinPeople])
# cs = ConstraintSet([XinPeople, YinPeople, X2inPeople, XeqY, X2eqY, X2eqX])
cs1 = ConstraintSet([X1inPeople, Y1inPeople])
cs_u = ConstraintSet([X2inPeople])

clause = ConstrainedClause([friendsXY, dislikesXY], [X2, X, Y], cs)
clause1 = ConstrainedClause([~friendsX1Y1, likesX1Y1], [X1, Y1], cs1)
unitclause = ConstrainedClause([friendsX2X2], [X2], cs_u)

cnf = CNF([clause, clause1, unitclause])
compiler = Compiler()
nnf = compiler.compile(cnf)

# print(unitclause.is_independent_from_other_clause(clause))


