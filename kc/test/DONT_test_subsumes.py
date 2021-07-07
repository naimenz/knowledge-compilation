"""Subsumption is hard.
 This file contains a number of test cases for figuring out how it should be calculated
 NOTE: I think this file is useless because I've removed is_subsumed_by_clause"""

from kc.data_structures import *

X, Y =  LogicalVariable('X'), LogicalVariable('Y')
Z, W = LogicalVariable('Z'), LogicalVariable('W')
X1, Y1 = LogicalVariable('X1'), LogicalVariable('Y1')
Z1, W1 =  LogicalVariable('Z1'), LogicalVariable('W1')
a, b, c, d = Constant('a'), Constant('b'), Constant('c'), Constant('d')

ab = SetOfConstants([a, b])
cd = SetOfConstants([c, d])

Xinab = InclusionConstraint(X, ab)
Yinab = InclusionConstraint(Y, ab)
Zinab = InclusionConstraint(Z, ab)

X1inab = InclusionConstraint(X1, ab)
Y1inab = InclusionConstraint(Y1, ab)
Z1inab = InclusionConstraint(Z1, ab)

Xincd = InclusionConstraint(X, cd)
Yincd = InclusionConstraint(Y, cd)

Xeqa = InclusionConstraint(X, SetOfConstants([a]))
Yeqa = InclusionConstraint(Y, SetOfConstants([a]))

XeqY = EqualityConstraint(X, Y)
XeqZ = EqualityConstraint(X, Z)
YeqZ = EqualityConstraint(Y, Z)
ZeqW = EqualityConstraint(Z, W)

p = Predicate('p', 2)
pXY = Literal(Atom(p, [X, Y]))
pX1Y1 = Literal(Atom(p, [X1, Y1]))
pXZ = Literal(Atom(p, [X, Z]))
literals = [pXY]

q = Predicate('q', 1)
qX = Literal(Atom(q, [X]))


clause1 = ConstrainedClause([pXY], [X, Y], ConstraintSet([Xinab, Yinab]))
clause1p = ConstrainedClause([pX1Y1], [X1, Y1], ConstraintSet([X1inab, Y1inab]))
clause2 = ConstrainedClause(literals, [X, Y], ConstraintSet([Xinab, Yinab, XeqY]))
clause3 = ConstrainedClause(literals, [X, Y], ConstraintSet([Xinab, Yinab, ~XeqY]))
def test_simple_subsumes():
    assert(clause1.is_subsumed_by_clause(clause1))
    assert(clause2.is_subsumed_by_clause(clause1))
    assert(not clause1.is_subsumed_by_clause(clause2))
    assert(clause3.is_subsumed_by_clause(clause1))
    assert(not clause1.is_subsumed_by_clause(clause3))

clause4 = ConstrainedClause(literals, [X, Y], ConstraintSet([Xinab, Yincd]))
clause5 = ConstrainedClause(literals, [X, Y], ConstraintSet([Xinab, Yincd, XeqY]))
clause6 = ConstrainedClause(literals, [X, Y], ConstraintSet([Xinab, Yincd, XeqZ]))
def test_disjoint_domains():
    assert(clause6.is_subsumed_by_clause(clause4))
    assert(not clause4.is_subsumed_by_clause(clause6))

fclause1 = ConstrainedClause([pXY], [X, Y], ConstraintSet([Xinab, Yinab]))
fclause2 = ConstrainedClause([pXZ], [X], ConstraintSet([Xinab]))
fclause3 = ConstrainedClause([pXY], [X, Y], ConstraintSet([Xinab, Yinab, YeqZ]))
fclause4 = ConstrainedClause([pXY], [X, Y], ConstraintSet([Xinab, Yinab, ~ZeqW]))
fclause4a = ConstrainedClause([pX1Y1], [X1, Y1], ConstraintSet([X1inab, Y1inab, ~ZeqW]))
def test_free_variables():
    assert(not fclause2.is_subsumed_by_clause(fclause1))
    assert(not fclause1.is_subsumed_by_clause(fclause2))
    assert(fclause3.is_subsumed_by_clause(fclause1))
    assert(not fclause2.is_subsumed_by_clause(fclause3))
    assert(fclause1.is_subsumed_by_clause(fclause1))
    assert(fclause4.is_subsumed_by_clause(fclause4a))

a_gamma = ConstrainedClause([qX], [X], ConstraintSet([Xinab]))
A = ConstrainedClause([qX], [X], ConstraintSet([Xinab, ~YeqZ]))
def test_example_42():
    # example 4.2 from splitting in the PhD
    assert(not a_gamma.is_subsumed_by_clause(A))
    assert(A.is_subsumed_by_clause(a_gamma))

def test_subsumes_self():
    assert(clause1.is_subsumed_by_clause(clause1))
    assert(clause1p.is_subsumed_by_clause(clause1))
    assert(clause1.is_subsumed_by_clause(clause1p))

    assert(clause2.is_subsumed_by_clause(clause2))
    assert(clause3.is_subsumed_by_clause(clause3))
    assert(clause4.is_subsumed_by_clause(clause4))
    # clause5 doesn't subsume itself because its cs unsatisfiable
    # i'm not sure how to feel about that
    # assert(clause5.is_subsumed_by_clause(clause5))
    assert(clause6.is_subsumed_by_clause(clause6))
    assert(fclause1.is_subsumed_by_clause(fclause1))
    assert(fclause2.is_subsumed_by_clause(fclause2))
    assert(fclause3.is_subsumed_by_clause(fclause3))
    assert(a_gamma.is_subsumed_by_clause(a_gamma))
    assert(A.is_subsumed_by_clause(A))
        
assert(fclause4.is_subsumed_by_clause(fclause4a))
