from kc.data_structures import *
from kc.compiler import *

X = LogicalVariable('X')
Y = LogicalVariable('Y')
X1 = LogicalVariable('X1')
Y1 = LogicalVariable('Y1')


alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
dislikes = Predicate('dislikes', 2)
likes = Predicate('likes', 2)
fun = Predicate('fun', 1)

friendsXY = Literal(Atom(friends, [X, Y]), True)
dislikesXY = Literal(Atom(dislikes, [X, Y]), True)
friendsX1Y1 = Literal(Atom(friends, [X1, Y1]), True)
funX1 = Literal(Atom(fun, [X1]), True)

People = SetOfConstants([alice, bob, charlie])
XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)
XeqY = EqualityConstraint(X, Y)
X1eqY1 = EqualityConstraint(X1, Y1)

cs1 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
cs2 = ConstraintSet([X1inPeople, Y1inPeople, ~X1eqY1])
clause1 = ConstrainedClause([dislikesXY, friendsXY], [X, Y], cs1)
clause2 = ConstrainedClause([funX1, ~friendsX1Y1], [X1, Y1], cs2)
cnf = CNF([clause1, clause2])
cnf.shattered = True # hack for now

compiler = Compiler()
nnf = compiler.compile(cnf)
