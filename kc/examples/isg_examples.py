from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf

X = LogicalVariable('X')
Y = LogicalVariable('Y')
X1 = LogicalVariable('X1')
Y1 = LogicalVariable('Y1')

Yfree = FreeVariable('Yfree')


alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
dislikes = Predicate('dislikes', 2)
likes = Predicate('likes', 2)
fun = Predicate('fun', 1)

friendsXY = Literal(Atom(friends, [X, Y]), True)
friendsXYfree = Literal(Atom(friends, [X, Yfree]), True)
friendsYX = Literal(Atom(friends, [Y, X]), True)
friendsYfreeX = Literal(Atom(friends, [Yfree, X]), True)
friendsYfreeYfree = Literal(Atom(friends, [Yfree, Yfree]), True)
friendsY1X1 = Literal(Atom(friends, [Y1, X1]), True)
friendsX1Y1 = Literal(Atom(friends, [X1, Y1]), True)
friendsYY = Literal(Atom(friends, [Y, Y]), True)

dislikesXY = Literal(Atom(dislikes, [X, Y]), True)
dislikesXYfree = Literal(Atom(dislikes, [X, Yfree]), True)
dislikesYfreeYfree = Literal(Atom(dislikes, [Yfree, Yfree]), True)
dislikesYY = Literal(Atom(dislikes, [Y, Y]), True)

funX = Literal(Atom(fun, [X]), True) 
funY = Literal(Atom(fun, [Y]), True) 
funY1 = Literal(Atom(fun, [Y1]), True) 
funX1 = Literal(Atom(fun, [X1]), True)

People = RootDomain([alice, bob, charlie], 'People')
XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)
XeqY = EqualityConstraint(X, Y)
XeqYfree = EqualityConstraint(X, Yfree)
X1eqY1 = EqualityConstraint(X1, Y1)

cs1 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
cs2 = ConstraintSet([X1inPeople, Y1inPeople, ~X1eqY1])
# cs2 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
clause1 = ConstrainedClause([dislikesXY, friendsXY], [X, Y], cs1)
# clause2 = ConstrainedClause([funY, ~friendsYX], [X, Y], cs2)
# clause2 = ConstrainedClause([funX1, ~friendsX1Y1], [X1, Y1], cs2)
clause2 = ConstrainedClause([funY1, ~friendsX1Y1], [X1, Y1], cs2)
cnf = CNF([clause1, clause2])
# cnf.shattered = True # hack for now

# c = ConstrainedClause([dislikesXY], [X, Y], ConstraintSet([XinPeople]))
# # print(c.substitute(Substitution([(Y, Yfree)])))
# gamma = ConstrainedClause([dislikesXYfree, friendsXYfree], [X], ConstraintSet([XinPeople]))
# # gamma = UnconstrainedClause([dislikesYfreeYfree, friendsYfreeYfree])
# # gamma = ConstrainedClause([dislikesXYfree, friendsXYfree], [X], ConstraintSet([XinPeople, ~XeqYfree]))
# print(gamma)
# A = ConstrainedAtom([friendsYfreeX], [X], ConstraintSet([XinPeople]))
# a_gamma = ConstrainedAtom([friendsYfreeYfree], [], ConstraintSet([]))
# print(A.subsumes(a_gamma))
# print(UnitPropagation.split(gamma, A))


compiler = Compiler()

nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)
