from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf
from kc.parsing import write_nnf_to_txt

X = LogicalVariable('X')
Y = LogicalVariable('Y')
X1 = LogicalVariable('X1')
Y1 = LogicalVariable('Y1')


alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
dislikes = Predicate('dislikes', 2)
fun = Predicate('fun', 1)

friendsXY = Literal(Atom(friends, [X, Y]))
dislikesXY = Literal(Atom(dislikes, [X, Y])) 

friendsX1Y1 = Literal(Atom(friends, [X1, Y1]))
funX1 = Literal(Atom(fun, [X1]))

People = RootDomain([alice, bob, charlie], 'People')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)

XeqY = EqualityConstraint(X, Y)
X1eqY1 = EqualityConstraint(X1, Y1)

cs1 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
cs2 = ConstraintSet([X1inPeople, Y1inPeople, ~X1eqY1])

clause1 = ConstrainedClause([friendsXY, dislikesXY], [X, Y], cs1)
clause2 = ConstrainedClause([funX1, ~friendsX1Y1], [X1, Y1], cs2)

cnf = CNF([clause1, clause2])
# cnf.shattered = True  # hack so I can test AC directly like in the PhD, when really should shatter first
compiler = Compiler()
nnf = compiler.compile(cnf)

current_node = nnf
# draw_nx_graph_from_nnf(nnf)
smoothed_nnf = nnf.do_smoothing(cnf)
# draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, 'TEST')
