"""Tests for the whole 'subsumes' set up for ConstrainedAtoms"""
from kc.data_structures import *

def test_is_not_trivial():
    X = LogicalVariable('X')
    Y = LogicalVariable('Y') alice = Constant('alice')
    bob = Constant('bob')
    charlie = Constant('charlie')
    dog = Constant('dog')
    cat = Constant('cat')
    People = SetOfConstants([alice, bob, charlie])
    Smokers = SetOfConstants([alice, bob])
    Animals = SetOfConstants([dog, cat])

    p = Predicate('p', 2)
    pXY = Literal(Atom(p, [X, Y]))

    XinPeople = InclusionConstraint(X, People)
    YinSmokers = InclusionConstraint(Y, Smokers)
    YinAnimals = InclusionConstraint(Y, Animals)
    XeqY = EqualityConstraint(X, Y)
    XneqY = ~XeqY
    cs_smokers = ConstraintSet([XinPeople, YinSmokers, XneqY])

    cs_animals = ConstraintSet([XinPeople, YinAnimals, XneqY])

    c_atom_smokers = ConstrainedAtom([pXY], [X, Y], cs_smokers)
    assert(XneqY.is_not_trivial(c_atom_smokers))

    c_atom_animals = ConstrainedAtom([pXY], [X, Y], cs_animals)
    assert(not XneqY.is_not_trivial(c_atom_animals))


def test_get_bound_variable_inequalities():
    X = LogicalVariable('X')
    Y = LogicalVariable('Y')
    Z = LogicalVariable('Z')
    alice = Constant('alice')
    bob = Constant('bob')
    charlie = Constant('charlie')
    dog = Constant('dog')
    cat = Constant('cat')
    People = SetOfConstants([alice, bob, charlie])
    Smokers = SetOfConstants([alice, bob])
    Animals = SetOfConstants([dog, cat])

    p = Predicate('p', 2)
    pXY = Literal(Atom(p, [X, Y]))
    pXZ = Literal(Atom(p, [X, Z]))

    XeqY = InequalityConstraint(X, Y)
    XeqZ = InequalityConstraint(X, Z)
    YeqZ = InequalityConstraint(Y, Z)

    cs = ConstraintSet([~XeqY, ~XeqZ, ~YeqZ])


    cclause1 = ConstrainedClause([pXY, pXZ], [X, Y], cs)
    target1 = set([~XeqY])

    cclause2 = ConstrainedClause([pXY, pXZ], [X, Z], cs)
    target2 = set([~XeqZ])

    cclause3 = ConstrainedClause([pXY, pXZ], [X, Y, Z], cs)
    target3 = set([~XeqY, ~XeqZ, ~YeqZ])

    assert(cclause1.get_bound_variable_inequalities() == target1)
    assert(cclause2.get_bound_variable_inequalities() == target2)
    assert(cclause3.get_bound_variable_inequalities() == target3)

print(test_get_bound_variable_inequalities())
