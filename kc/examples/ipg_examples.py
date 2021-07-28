"""This is the example from the IPG section, Example 4.15"""
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf

X = LogicalVariable('X')
Y = LogicalVariable('Y')
X1 = LogicalVariable('X1')
Y1 = LogicalVariable('Y1')


alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
enemies = Predicate('enemies', 2)

friendsXY = Literal(Atom(friends, [X, Y]), True)
friendsY1X1 = Literal(Atom(friends, [Y1, X1]), True)
enemiesXY = Literal(Atom(enemies, [X, Y]), True)
enemiesX1Y1 = Literal(Atom(enemies, [X1, Y1]), True)

# People = SetOfConstants([alice, bob, charlie])
People = RootDomain([alice, bob, charlie], 'People')
XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)

cs1 = ConstraintSet([XinPeople, YinPeople])
cs2 = ConstraintSet([X1inPeople, Y1inPeople])
clause1 = ConstrainedClause([~friendsXY, ~enemiesXY], [X, Y], cs1)
clause2 = ConstrainedClause([~friendsY1X1, ~enemiesX1Y1], [X1, Y1], cs2)
cnf = CNF([clause1, clause2])

compiler = Compiler()
# cnf.shattered = True  # hack to make it more like PhD
nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.get_smoothed_node(cnf)
# smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)
