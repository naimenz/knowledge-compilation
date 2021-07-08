"""Tests for the whole 'subsumes' set up for ConstrainedAtoms"""
from kc.data_structures import *

def test_is_not_trivial():
    X = LogicalVariable('X')
    Y = LogicalVariable('Y')
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
    assert (not XneqY.is_not_trivial(c_atom_animals))

