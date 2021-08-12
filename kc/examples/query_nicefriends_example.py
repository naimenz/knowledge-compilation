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
f1 = Predicate('f_1', 2)

friendsXY = Literal(Atom(friends, [X, Y]))
niceY = Literal(Atom(nice, [Y]))
f1XY = Literal(Atom(f1, [X, Y]))
niceguy = Literal(Atom(nice, [guy]))

People = RootDomain([alice, bob, charlie, guy], 'person')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([~niceY, friendsXY], [X, Y], cs)
auxiliary_clauses = make_auxiliary_predicates_for_clauses([clause])
query = UnconstrainedClause([niceguy])

cnf = CNF(auxiliary_clauses + [query])
print(cnf)
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()

nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, 'nice_query')
