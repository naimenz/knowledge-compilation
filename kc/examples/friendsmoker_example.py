"""This is the friendsmoker example from Ex 3.17
In non-CNF form, this is smokes(X) ^ friends(X, Y) => smokes(Y)"""
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf

X = LogicalVariable('X')
Y = LogicalVariable('Y')

X1 = LogicalVariable('X1')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
smokes = Predicate('smokes', 1)

friendsXX = Literal(Atom(friends, [X, X]))
friendsXY = Literal(Atom(friends, [X, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))

friendsXX1 = Literal(Atom(friends, [X, X1]))
smokesX = Literal(Atom(smokes, [X]))
smokesX1 = Literal(Atom(smokes, [X1]))

People = RootDomain([alice, bob, charlie], 'People')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
XeqY = EqualityConstraint(X, Y)

cs = ConstraintSet([XinPeople, YinPeople])
# cs = ConstraintSet([XinPeople, YinPeople, ~XeqY])
csX = ConstraintSet([XinPeople])

clause = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs)
clauseXX = ConstrainedClause([~smokesX, ~friendsXX, smokesX], [X], csX)

D = DomainVariable('D', People)
Dcomp = D.complement

X1inPeople = InclusionConstraint(X1, People)

XinD = InclusionConstraint(X, D)
XinDcomp = InclusionConstraint(X, Dcomp)

X1inD = InclusionConstraint(X1, D)
X1inDcomp = InclusionConstraint(X1, Dcomp)

gamma_s1 = ConstrainedClause([~smokesX, ~friendsXX1, smokesX1], [X, X1],
        ConstraintSet([XinPeople, X1inPeople, XinD, X1inDcomp]))

gamma_s2 = ConstrainedClause([~smokesX, ~friendsXX1, smokesX1], [X, X1],
        ConstraintSet([XinPeople, X1inPeople, X1inD, XinDcomp]))

c_literal = UnitClause([smokesX1], [X, X1],
        ConstraintSet([XinPeople, X1inPeople, X1inD, XinDcomp]))

unit_clause = ConstrainedAtom([smokesX], [X],
        ConstraintSet([XinPeople, XinD]))
# print(c_literal.is_subsumed_by_literal(unit_clause))
# print(UnitPropagation.condition(gamma_s2, unit_clause))

# clause = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs.join(ConstraintSet([EqualityConstraint(X, Y)])))


cnf = CNF([clause])
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()
nnf = compiler.compile(cnf)
smoothed_nnf = nnf.get_smoothed_node()

draw_nx_graph_from_nnf(nnf)

draw_nx_graph_from_nnf(smoothed_nnf)
