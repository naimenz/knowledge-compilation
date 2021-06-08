from kc.data_structures.logicalterms import *

def test_constant_equality():
    a1 = Constant('a')
    a2 = Constant('a')
    b  = Constant('b')


    assert(a1 == a2)
    assert(a1 != b)
    assert(b == b)

def test_variable_equality():
    X1 = LogicalVariable('X')
    X2 = LogicalVariable('X')
    Y  = LogicalVariable('Y')

    assert(X1 == X2)
    assert(X1 != Y)
    assert(Y == Y)

def test_cross_equality():
    a = Constant('a')
    var_a = LogicalVariable('a')
    assert(var_a != a)


