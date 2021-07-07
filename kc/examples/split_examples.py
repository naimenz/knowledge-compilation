from kc.data_structures import *
from kc.compiler import *

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

p = Predicate('p', 1)
q = Predicate('q', 1)

pX = Literal(Atom(p, [X]))
qX = Literal(Atom(q, [X]))
pX1 = Literal(Atom(p, [X1]))

People = SetOfConstants([alice, bob, charlie])

XinPeople = InclusionConstraint(X, People)
X1inPeople = InclusionConstraint(X1, People)

YeqZ = EqualityConstraint(Y, Z)

cs = ConstraintSet([XinPeople])
cs1 = ConstraintSet([X1inPeople, ~YeqZ])

clause = ConstrainedClause([pX, qX], [X], cs)
c_atom = ConstrainedAtom([pX1], [X1], cs1)

print(UnitPropagation.split(clause, c_atom))

# print(unitclause.is_independent_from_other_clause(clause))


