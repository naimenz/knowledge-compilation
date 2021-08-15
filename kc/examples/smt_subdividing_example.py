
"""This is the friendsmoker example from Ex 3.17
In non-CNF form, this is f1(X, Y) <=> ( smokes(X) ^ friends(X, Y) => smokes(Y) )"""
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf
from kc.parsing import write_nnf_to_txt, make_auxiliary_predicates_for_clauses

X = LogicalVariable('X')
Y = LogicalVariable('Y')

alice = Constant('Alice')
bob = Constant('Bob')
charlie = Constant('Charlie')
guy = Constant('Guy')

friends = Predicate('friends', 2)
smokes = Predicate('smokes', 1)
age_10_20 = SMTPredicate('age', 1, 10, 20)
age_15_30 = SMTPredicate('age', 1, 15, 30)

friendsXY = Literal(Atom(friends, [X, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))
age_10_20X = Literal(Atom(age_10_20, [X]))
age_15_30Y = Literal(Atom(age_15_30, [Y]))

People = RootDomain([alice, bob, charlie, guy], 'person')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause1 = ConstrainedClause([~smokesX, ~friendsXY, age_10_20X], [X, Y], cs)
clause2 = ConstrainedClause([~smokesX, ~friendsXY, age_15_30Y], [X, Y], cs)

cnf = CNF([clause1, clause2])
print(cnf)
print(cnf.subdivide_ranges())
