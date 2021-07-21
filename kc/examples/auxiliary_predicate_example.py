"""Testing the new functions for auxiliary predicate generation"""
from kc.data_structures import *
from kc.parsing import make_auxiliary_predicate_for_clause

X = LogicalVariable('X')
Y = LogicalVariable('Y')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
smokes = Predicate('smokes', 1)
f1 = Predicate('f1', 2)

friendsXY = Literal(Atom(friends, [X, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))
f1XY = Literal(Atom(f1, [X, Y]))

People = RootDomain([alice, bob, charlie], 'People')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([smokesY, ~smokesX, ~friendsXY], [X, Y], cs)

clause1 = ConstrainedClause([~f1XY, smokesY, ~smokesX, ~friendsXY], [X, Y], cs)
clause2 = ConstrainedClause([f1XY, ~smokesY], [X, Y], cs)
clause3 = ConstrainedClause([f1XY, friendsXY], [X, Y], cs)
clause4 = ConstrainedClause([f1XY, smokesX], [X, Y], cs)

cnf = CNF([clause1, clause2, clause3, clause4])

print(f'Hand-made CNF:\n{cnf}')
print(f'Generated CNF:\n{CNF(make_auxiliary_predicate_for_clause(clause))}')

