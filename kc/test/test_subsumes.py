"""Tests for the whole 'subsumes' set up for ConstrainedAtoms
NOTE: The does_not_subsume tests only work when the function is 
edited to return debug strings of '1', '2', '3', '4'."""
from kc.data_structures import *
X = LogicalVariable('X')
Y = LogicalVariable('Y')
Z = LogicalVariable('Z')
W = LogicalVariable('W')

F = LogicalVariable('F')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

dog = Constant('dog')
cat = Constant('cat')
People = SetOfConstants([alice, bob, charlie])
Smokers = SetOfConstants([alice, bob])
Animals = SetOfConstants([dog, cat])

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
ZinPeople = InclusionConstraint(Z, People)
WinPeople = InclusionConstraint(W, People)

YinAnimals = InclusionConstraint(Y, Animals)
WinAnimals = InclusionConstraint(W, Animals)

p = Predicate('p', 2)
pXY = Literal(Atom(p, [X, Y]))
pZW = Literal(Atom(p, [Z, W]))
pXZ = Literal(Atom(p, [X, Z]))
pXX = Literal(Atom(p, [X, X]))
pYY = Literal(Atom(p, [Y, Y]))
pZZ = Literal(Atom(p, [Z, Z]))
pXa = Literal(Atom(p, [X, alice]))
pZa = Literal(Atom(p, [Z, alice]))

XeqY = EqualityConstraint(X, Y)
XeqZ = EqualityConstraint(X, Z)
YeqZ = EqualityConstraint(Y, Z)

XeqF = EqualityConstraint(X, F)
ZeqF = EqualityConstraint(Z, F)

Xeqa = InclusionConstraint(X, SetOfConstants([alice]))
Yeqa = InclusionConstraint(Y, SetOfConstants([alice]))
Zeqa = InclusionConstraint(Z, SetOfConstants([alice]))
Weqa = InclusionConstraint(W, SetOfConstants([alice]))

def test_is_not_trivial():

    XinPeople = InclusionConstraint(X, People)
    YinSmokers = InclusionConstraint(Y, Smokers)
    YinAnimals = InclusionConstraint(Y, Animals)
    XneqY = ~XeqY
    cs_smokers = ConstraintSet([XinPeople, YinSmokers, XneqY])

    cs_animals = ConstraintSet([XinPeople, YinAnimals, XneqY])

    c_atom_smokers = ConstrainedAtom([pXY], [X, Y], cs_smokers)
    assert(XneqY.is_not_trivial(c_atom_smokers))

    c_atom_animals = ConstrainedAtom([pXY], [X, Y], cs_animals)
    assert(not XneqY.is_not_trivial(c_atom_animals))


def test_get_bound_variable_inequalities():
    cs = ConstraintSet([~XeqY, ~XeqZ, ~YeqZ])


    cclause1 = ConstrainedClause([pXY, pXZ], [X, Y], cs)
    target1 = set([~XeqY])

    cclause2 = ConstrainedClause([pXY, pXZ], [X, Z], cs)
    target2 = set([~XeqZ])

    cclause3 = ConstrainedClause([pXY, pXZ], [X, Y, Z], cs)
    target3 = set([~XeqY, ~XeqZ, ~YeqZ])

    cclause4 = ConstrainedClause([pXY, pXZ], [Z, W], ConstraintSet([XeqY]))
    target4 = set()

    assert(cclause1.get_bound_variable_inequalities() == target1)
    assert(cclause2.get_bound_variable_inequalities() == target2)
    assert(cclause3.get_bound_variable_inequalities() == target3)
    assert(cclause4.get_bound_variable_inequalities() == target4)

def test_get_constant_or_free_inequalities():
    X = LogicalVariable('X')
    Y = LogicalVariable('Y')
    Z = LogicalVariable('Z')
    W = LogicalVariable('W')
    alice = Constant('alice')
    bob = Constant('bob')
    charlie = Constant('charlie')

    People = SetOfConstants([alice, bob, charlie])

    p = Predicate('p', 2)
    pXY = Literal(Atom(p, [X, Y]))
    pXZ = Literal(Atom(p, [X, Z]))

    Xeqa = InclusionConstraint(X, SetOfConstants([alice]))
    Yeqb = InclusionConstraint(Y, SetOfConstants([bob]))
    Zeqc = InclusionConstraint(Z, SetOfConstants([charlie]))
    XinPeople = InclusionConstraint(X, People)

    cs = ConstraintSet([~Xeqa, ~Yeqb, ~Zeqc, ~XinPeople])

    cclause1 = ConstrainedClause([pXY, pXZ], [X, Y], cs)
    target1 = set([~Xeqa, ~Yeqb])

    cclause2 = ConstrainedClause([pXY, pXZ], [X, Z], cs)
    target2 = set([~Xeqa, ~Zeqc])

    cclause3 = ConstrainedClause([pXY, pXZ], [X, Y, Z], cs)
    target3 = set([~Xeqa, ~Yeqb, ~Zeqc])

    cclause4 = ConstrainedClause([pXY, pXZ], [W], cs)
    target4 = set()

    assert(cclause1.get_constant_or_free_inequalities() == target1)
    assert(cclause2.get_constant_or_free_inequalities() == target2)
    assert(cclause3.get_constant_or_free_inequalities() == target3)
    assert(cclause4.get_constant_or_free_inequalities() == target4)

def test_does_not_subsume1():
    catom1 = ConstrainedAtom([pXa], [X], ConstraintSet([XinPeople]))
    catom2 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople]))

    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom2)
    assert(catom1.does_not_subsume(catom2, mgu) == "1")
    assert(catom2.does_not_subsume(catom1, mgu) != "1")

    catom3 = ConstrainedAtom([pXY], [X, Y], ConstraintSet([XinPeople, YinPeople]))
    mgu = catom2.get_constrained_atom_mgu_eq_classes(catom3)
    assert(catom2.does_not_subsume(catom3, mgu) != "1")

def test_does_not_subsume1_free_vars():
    catom1 = ConstrainedAtom([pXY], [X], ConstraintSet([XinPeople]))
    catom2 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople]))

    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom2)
    assert(catom1.does_not_subsume(catom2, mgu) == "1")
    assert(catom2.does_not_subsume(catom1, mgu) != "1")

    catom3 = ConstrainedAtom([pZa], [Z], ConstraintSet([ZinPeople]))

    # NOTE: I'm not sure what should happen here but it passes
    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom3)
    assert(catom3.does_not_subsume(catom1, mgu) != "1")
    assert(catom1.does_not_subsume(catom3, mgu) != "1")


def test_does_not_subsume2():
    catom1 = ConstrainedAtom([pXX], [X], ConstraintSet([XinPeople]))
    catom2 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople]))

    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom2)
    assert(catom1.does_not_subsume(catom2, mgu) == "2")
    assert(catom2.does_not_subsume(catom1, mgu) != "2")

    catom3 = ConstrainedAtom([pYY], [Y], ConstraintSet([XinPeople, YinPeople]))
    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom3)
    assert(catom1.does_not_subsume(catom3, mgu) != "2")

def test_does_not_subsume3():
    catom1 = ConstrainedAtom([pXY], [X, Y], ConstraintSet([XinPeople, YinPeople, ~Xeqa]))
    catom2 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople]))

    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom2)
    assert(catom1.does_not_subsume(catom2, mgu) == "3")
    assert(catom2.does_not_subsume(catom1, mgu) != "3")

    catom3 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople, ~Zeqa]))
    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom3)
    assert(catom1.does_not_subsume(catom3, mgu) != "3")
    assert(catom3.does_not_subsume(catom1, mgu) != "3")

def test_does_not_subsume3_free_vars():
    catom1 = ConstrainedAtom([pXY], [X, Y], ConstraintSet([XinPeople, YinPeople, ~XeqF]))
    catom2 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople]))

    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom2)
    assert(catom1.does_not_subsume(catom2, mgu) == "3")
    assert(catom2.does_not_subsume(catom1, mgu) != "3")

    catom3 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople, ~ZeqF]))
    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom3)
    assert(catom1.does_not_subsume(catom3, mgu) != "3")
    assert(catom3.does_not_subsume(catom1, mgu) != "3")

def test_does_not_subsume4():
    catom1 = ConstrainedAtom([pXY], [X, Y], ConstraintSet([XinPeople, YinPeople, ~XeqY]))
    catom2 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinPeople]))

    mgu = catom1.get_constrained_atom_mgu_eq_classes(catom2)
    assert(catom1.does_not_subsume(catom2, mgu) == "4")
    assert(catom2.does_not_subsume(catom1, mgu) != "4")

    catom3 = ConstrainedAtom([pXY], [X, Y], ConstraintSet([XinPeople, YinAnimals, ~XeqY]))
    catom4 = ConstrainedAtom([pZW], [Z, W], ConstraintSet([ZinPeople, WinAnimals]))
    mgu = catom3.get_constrained_atom_mgu_eq_classes(catom4)
    assert(catom3.does_not_subsume(catom4, mgu) != "4")
