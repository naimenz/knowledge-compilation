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
f1 = Predicate('f_1', 2)

friendsXY = Literal(Atom(friends, [X, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))
f1XY = Literal(Atom(f1, [X, Y]))

People = RootDomain([alice, bob, charlie, guy], 'person')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs)
# clause1 = ConstrainedClause([~f1XY, smokesY, ~smokesX, ~friendsXY], [X, Y], cs)
# clause2 = ConstrainedClause([f1XY, ~smokesY], [X, Y], cs)
# clause3 = ConstrainedClause([f1XY, friendsXY], [X, Y], cs)
# clause4 = ConstrainedClause([f1XY, smokesX], [X, Y], cs)
auxiliary_clauses = make_auxiliary_predicates_for_clauses([clause])

cnf = CNF(auxiliary_clauses)
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()

import time
tic = time.perf_counter()
N = 10
for i in range(N):
    nnf = compiler.compile(cnf)
    # draw_nx_graph_from_nnf(nnf)

    smoothed_nnf = nnf.do_smoothing(cnf)
toc = time.perf_counter()
print(f'10 iterations took {toc - tic:0.4f} seconds')
draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, 'auxiliary_theory')
