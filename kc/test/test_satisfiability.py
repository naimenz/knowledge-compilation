"""Satisfiability is hard.
 This file contains a number of test cases for figuring out how it should be calculated"""

from kc.data_structures import *

X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
a, b, c, d = Constant('a'), Constant('b'), Constant('c'), Constant('d')

ab = SetOfConstants([a, b])
cd = SetOfConstants([c, d])

Xinab = InclusionConstraint(X, ab)
Yinab = InclusionConstraint(Y, ab)
Zinab = InclusionConstraint(Z, ab)

Xincd = InclusionConstraint(X, cd)
Yincd = InclusionConstraint(Y, cd)

Xeqa = InclusionConstraint(X, SetOfConstants([a]))
Yeqa = InclusionConstraint(Y, SetOfConstants([a]))

XeqY = EqualityConstraint(X, Y)
XeqZ = EqualityConstraint(X, Z)
YeqZ = EqualityConstraint(Y, Z)

def test_bound_variables():
    assert(ConstraintSet([Xinab, Yinab]).is_satisfiable())
    assert(ConstraintSet([Xinab, Yinab, XeqY]).is_satisfiable())
    assert(ConstraintSet([Xinab, Yinab, ~XeqY]).is_satisfiable())

    assert(ConstraintSet([Xinab, Yincd]).is_satisfiable())
    assert(not ConstraintSet([Xinab, Yincd, XeqY]).is_satisfiable())
    assert(ConstraintSet([Xinab, Yincd, ~XeqY]).is_satisfiable())

def test_free_variables():
    assert(ConstraintSet([Xinab, Yincd, XeqZ]).is_satisfiable())
    assert(not ConstraintSet([Xinab, Yincd, XeqZ, YeqZ]).is_satisfiable())

def test_constant_equality():
    assert(ConstraintSet([Xeqa, Yeqa, XeqY]).is_satisfiable())
    assert(not ConstraintSet([Xeqa, Yeqa, ~XeqY]).is_satisfiable())

def test_mutual_inequality():
    assert(not ConstraintSet([Xinab, Yinab, Zinab, ~XeqY, ~XeqZ, ~YeqZ]).is_satisfiable())

def test_false_constraint():
    assert(ConstraintSet([Xinab, Yinab, FalseConstraint("false constraint")]).is_satisfiable())
