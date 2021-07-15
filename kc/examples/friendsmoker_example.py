"""This is the friendsmoker example from Ex 3.17
In non-CNF form, this is smokes(X) ^ friends(X, Y) => smokes(Y)"""
from kc.data_structures import *
from kc.compiler import *

X = LogicalVariable('X')
Y = LogicalVariable('Y')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
smokes = Predicate('smokes', 1)

friendsXY = Literal(Atom(friends, [X, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))

People = RootDomain([alice, bob, charlie], 'People')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs)

cnf = CNF([clause])
cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
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
