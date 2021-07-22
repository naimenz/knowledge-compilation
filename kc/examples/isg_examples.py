from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf

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
funX = Literal(Atom(fun, [X]), True) 
friendsX1Y1 = Literal(Atom(friends, [X1, Y1]), True)
funX1 = Literal(Atom(fun, [X1]), True)

People = RootDomain([alice, bob, charlie], 'People')
XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)
XeqY = EqualityConstraint(X, Y)
X1eqY1 = EqualityConstraint(X1, Y1)

cs1 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
# cs2 = ConstraintSet([X1inPeople, Y1inPeople, ~X1eqY1])
cs2 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
clause1 = ConstrainedClause([dislikesXY, friendsXY], [X, Y], cs1)
# clause2 = ConstrainedClause([funX1, ~friendsX1Y1], [X1, Y1], cs2)
clause2 = ConstrainedClause([funX, ~friendsXY], [X, Y], cs2)
cnf = CNF([clause1, clause2])
cnf.shattered = True # hack for now

compiler = Compiler()
nnf = compiler.compile(cnf)

current_node = nnf
graph = build_nx_graph_from_nnf(current_node)
print(graph)

import matplotlib.pyplot as plt
import pydot
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
pos = graphviz_layout(graph, prog="dot")
nx.draw(graph, pos)

label_pos = {}
y_off = 10  # offset on the y axis

for k, v in pos.items():
    label_pos[k] = (v[0], v[1]+y_off)
node_labels = nx.get_node_attributes(graph, 'label')
nx.draw_networkx_labels(graph, label_pos, labels=node_labels)
plt.show()
