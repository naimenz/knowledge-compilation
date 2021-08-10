"""This is the friendsmoker example from Ex 3.17
In non-CNF form, this is f1(X, Y) <=> ( smokes(X) ^ friends(X, Y) => smokes(Y) )"""
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf
from kc.parsing import write_nnf_to_txt

X = LogicalVariable('X')
Y = LogicalVariable('Y')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')
guy = Constant('guy')

friends = Predicate('friends', 2)
smokes = Predicate('smokes', 1)
f1 = Predicate('f_1', 2)

friendsXY = Literal(Atom(friends, [X, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))
f1XY = Literal(Atom(f1, [X, Y]))

People = RootDomain([alice, bob, charlie, guy], 'person')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause1 = ConstrainedClause([~f1XY, smokesY, ~smokesX, ~friendsXY], [X, Y], cs)
clause2 = ConstrainedClause([f1XY, ~smokesY], [X, Y], cs)
clause3 = ConstrainedClause([f1XY, friendsXY], [X, Y], cs)
clause4 = ConstrainedClause([f1XY, smokesX], [X, Y], cs)

cnf = CNF([clause1, clause2, clause3, clause4])
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()

X1, Y1 = LogicalVariable('X1'), LogicalVariable('Y1')
friendsX1X = Literal(Atom(friends, [X1, X]))
X1inPeople = InclusionConstraint(X1, People)

# nc = ConstrainedClause([friendsX1X], [X1, X], ConstraintSet([X1inPeople, XinPeople]))
# nc1 = nc.rename_bound_variables(("X1", "Y1"))
# print(f'{nc = }')
# print(f'{nc1 = }')


nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, 'auxiliary_theory')
