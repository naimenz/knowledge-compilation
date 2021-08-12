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
nice = Predicate('nice', 1)
a1 = Predicate('a1', 2)

friendsXY = Literal(Atom(friends, [X, Y]))
friendsXX = Literal(Atom(friends, [X, X]))
niceX = Literal(Atom(nice, [X]))
niceY = Literal(Atom(nice, [Y]))
a1XY = Literal(Atom(a1, [X, Y]))
a1XX = Literal(Atom(a1, [X, X]))
niceguy = Literal(Atom(nice, [guy]))

People = RootDomain([alice, bob, charlie, guy], 'person')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
XeqY = EqualityConstraint(X, Y)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([~niceY, friendsXY], [X, Y], cs)
auxiliary_clauses = make_auxiliary_predicates_for_clauses([clause])
query = UnconstrainedClause([niceguy])

cs1 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
cs2 = ConstraintSet([XinPeople])

# # shattering independent test
# clause1 = ConstrainedClause([~a1XY, ~niceY, friendsXY], [X, Y], cs1)
# clause2 = ConstrainedClause([~a1XX, ~niceX, friendsXX], [X], cs2)
# clause3 = ConstrainedClause([~a1XY, friendsXY], [X, Y], cs1)
# clause4 = ConstrainedClause([~a1XX, friendsXX], [X], cs2)
# clause5 = ConstrainedClause([~a1XY, niceY], [X, Y], cs1)
# clause6 = ConstrainedClause([~a1XX, niceX], [X], cs2)
# print(clause5.is_independent_from_other_clause(clause4))


cnf = CNF(auxiliary_clauses)
print(cnf)
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()

nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, 'nice_theory')
