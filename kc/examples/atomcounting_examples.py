from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf

X = LogicalVariable('X')
Y = LogicalVariable('Y')
Z = LogicalVariable('Z')
X1 = LogicalVariable('X1')
Y1 = LogicalVariable('Y1')
Z1 = LogicalVariable('Z1')
X2 = LogicalVariable('X2')


alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')
diana = Constant('diana')

friends = Predicate('friends', 2)
fun = Predicate('fun', 1)

friendsXY = Literal(Atom(friends, [X, Y]))
friendsYX = Literal(Atom(friends, [Y, X]))
funX = Literal(Atom(fun, [X]))

funalice = Literal(Atom(fun, [alice]))
friendsaliceY = Literal(Atom(friends, [alice, Y]))

friendsY1X1 = Literal(Atom(friends, [Y1, X1]))
funX1 = Literal(Atom(fun, [X1]))

People = RootDomain([alice, bob, charlie], 'People')

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)

Xeqalice = InclusionConstraint(X, SetOfConstants([alice]))

D = DomainVariable('D', People)
# print(D.complement)
XinD = InclusionConstraint(X, D)
XinDcomp = InclusionConstraint(X, D.complement)
X1inD = InclusionConstraint(X1, D)
X1inDcomp = InclusionConstraint(X1, D.complement)

cs1 = ConstraintSet([XinPeople, YinPeople])
# cs1 = ConstraintSet([XinPeople, YinPeople, Xeqalice])
cs2 = ConstraintSet([X1inPeople, Y1inPeople])

# clause1 = ConstrainedClause([funalice, ~friendsaliceY], [X, Y], cs1)
clause1 = ConstrainedClause([funX, ~friendsXY], [X, Y], cs1)
# print(clause1)
# sub = Substitution([(X, alice)])
# print(sub)
# print(clause1.substitute(sub))
# print(clause1.propagate_equality_constraints())
clause2 = ConstrainedClause([funX1, ~friendsY1X1], [X1, Y1], cs2)
# clause2 = ConstrainedClause([funX, ~friendsYX], [X, Y], cs1)
c_atom = ConstrainedAtom([friendsXY], [X, Y], cs1)
other_c_atom = ConstrainedAtom([friendsYX], [X, Y], cs1)

d_other_c_atom = UnitPropagation._make_variables_different(c_atom, other_c_atom)
new_other_atom, new_clause = UnitPropagation._align_variables(c_atom, other_c_atom, clause2)
# print(f'atoms: {c_atom = }, {other_c_atom = }')
# print(f'{d_other_c_atom = }')
# print(f'{new_other_atom = }, {new_clause = }')

splitclause1 = ConstrainedClause([funX, ~friendsXY], [X, Y], cs1.join(ConstraintSet([XinD])))
splitclause2 = ConstrainedClause([funX1, ~friendsY1X1], [X1, Y1], cs2.join(ConstraintSet([X1inDcomp])))
u1 = ConstrainedAtom([funX], [X], ConstraintSet([XinPeople, XinD]))
u2 = UnitClause([~funX1], [X1], ConstraintSet([X1inPeople, X1inDcomp]))
bigger_uclause = ConstrainedAtom([funX], [X], ConstraintSet([XinPeople]))
# print(f'{splitclause1.is_subsumed_by_literal(u1) = }')
# print(f'{splitclause2.is_subsumed_by_literal(u2) = }')
# print(f'{ UnitPropagation.condition(splitclause1, u1) = }')
# print(f'{ UnitPropagation.condition(splitclause2, u2) = }')
# print(f'{u1=}')
# for ca in clause2.get_constrained_atoms():
#     print(f"================= {ca} =============")
#     print(ca.needs_splitting(u1))
# print(*UnitPropagation.split(clause2, u1), sep='\n')
cnf = CNF([clause1, clause2])
# cnf.shattered = True  # hack so I can test AC directly like in the PhD, when really should shatter first

compiler = Compiler()
nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.get_smoothed_node()
# smoothed_nnf = nnf.do_smoothing(cnf)
draw_nx_graph_from_nnf(smoothed_nnf)
