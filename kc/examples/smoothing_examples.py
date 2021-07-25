
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf

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
friendsX1Y1 = Literal(Atom(friends, [X1, Y1]))
funX = Literal(Atom(fun, [X]))

funalice = Literal(Atom(fun, [alice]))
friendsaliceY = Literal(Atom(friends, [alice, Y]))
friendsdianaY = Literal(Atom(friends, [diana, Y]))

friendsY1X1 = Literal(Atom(friends, [Y1, X1]))
funX1 = Literal(Atom(fun, [X1]))

People = RootDomain([alice, bob, charlie], 'People')
D = DomainVariable('D', People)

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

XinD = InclusionConstraint(X, D)
YinD = InclusionConstraint(Y, D)

X1inPeople = InclusionConstraint(X1, People)
Y1inPeople = InclusionConstraint(Y1, People)

XeqY = EqualityConstraint(X, Y)
Xeqalice = InclusionConstraint(X, SetOfConstants([alice]))

u1 = ConstrainedAtom([friendsXY], [X, Y], ConstraintSet([XinPeople, YinPeople]))
u2 = ConstrainedAtom([friendsXY], [X, Y], ConstraintSet([XinD, YinD]))
u3 = ConstrainedAtom([friendsdianaY], [Y], ConstraintSet([YinPeople]))
# print(u3.is_independent_from_other_clause(u2))
nnf = TrueNode()
# print(CNF(nnf._make_independent(set([u1, u2, u3]))))

c1 = ConstrainedAtom([funX], [X], ConstraintSet([XinPeople]))
c2 = ConstrainedAtom([funX], [X], ConstraintSet([XinD]))

print(nnf.A_without_B(set([u1, u3]), set([u1, u2, c1])))
print(nnf.A_without_B(set([u1, u3]), set([u1, c2])))
print(nnf.A_without_B(set([u1]), set([u1])))

