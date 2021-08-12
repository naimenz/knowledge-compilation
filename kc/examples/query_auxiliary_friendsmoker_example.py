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
friendsYX = Literal(Atom(friends, [Y, X]))
friendsGuyX = Literal(Atom(friends, [guy, X]))
friendsGuyY = Literal(Atom(friends, [guy, Y]))
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))
f1XY = Literal(Atom(f1, [X, Y]))
smokesguy = Literal(Atom(smokes, [guy]))

People = RootDomain([alice, bob, charlie, guy], 'person')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
XeqGuy = InclusionConstraint(X, SetOfConstants([guy]))
YeqGuy = InclusionConstraint(Y, SetOfConstants([guy]))
XeqY = EqualityConstraint(X, Y)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs)
auxiliary_clauses = make_auxiliary_predicates_for_clauses([clause])
query = UnconstrainedClause([smokesguy])

# possible = [ConstrainedAtom([friendsYX], [X, Y], cs)]
# circuit = [
#         ConstrainedAtom([friendsYX], [X, Y], ConstraintSet([XinPeople, YinPeople, ~XeqGuy, ~YeqGuy, ~XeqY])),
#         ConstrainedAtom([friendsGuyX], [X], ConstraintSet([XinPeople])),
#             ]
# possible = [ConstrainedAtom([friendsXY], [X, Y], cs)]
# circuit = [
#         ConstrainedAtom([friendsXY], [X, Y], ConstraintSet([XinPeople, YinPeople, ~XeqGuy, ~YeqGuy, ~XeqY])),
#         # ConstrainedAtom([friendsGuyX], [X], ConstraintSet([XinPeople])),
#             ]
# nnf = TrueNode()
# print(nnf.A_without_B(possible, circuit))
# possible = ConstrainedAtom([friendsXY], [X, Y], cs)
# circuit = ConstrainedAtom([friendsXY], [X, Y], ConstraintSet([XinPeople, YinPeople, ~XeqGuy, ~YeqGuy, ~XeqY]))
# print(f'{circuit.subtract_from_c_atom(possible) = }')


cnf = CNF(auxiliary_clauses + [query])
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()

nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)

# DEBUG trying unsmoothed
write_nnf_to_txt(smoothed_nnf, 'auxiliary_query')
