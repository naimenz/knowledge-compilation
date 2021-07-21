from kc.data_structures import *

likes = Predicate('likes', 2)
friends = Predicate('friends', 2)
X = LogicalVariable('X')
Y = LogicalVariable('Y')
a = Constant('a')
b = Constant('b')
c = Constant('c')
atom = Atom(likes, [X, a])
notlikesXa = Literal(Atom(likes, [X, a]), False)
friendsYY = Literal(Atom(friends, [Y, Y]), True)
uclause = UnconstrainedClause([notlikesXa, friendsYY])

D = SetOfConstants([a,b,c])
XeqY = EqualityConstraint(X, Y)
Xeqa = InclusionConstraint(X, SetOfConstants([a]))
Yeqa = InclusionConstraint(Y, SetOfConstants([a]))
Yneqa = NotInclusionConstraint(Y, SetOfConstants([a]))
XinD = InclusionConstraint(X, D)
YinD = InclusionConstraint(Y, D)

cs = ConstraintSet([Xeqa, Yneqa, XinD, YinD])
cclause = ConstrainedClause([notlikesXa, friendsYY], [X, Y], cs)
cnf = CNF([uclause, cclause])

subXa = Substitution([(X, a)])
subaX = Substitution([(a, X)])
subYX = Substitution([(Y, X)])

def test_atom():
    returned_atom = atom.substitute(subXa)
    target_atom = Atom(likes, [a, a])
    assert(returned_atom == target_atom)

def test_literal():
    returned_literal = notlikesXa.substitute(subXa)
    target_literal = Literal(Atom(likes, [a, a]), False)
    assert(returned_literal == target_literal)

def test_uclause():
    returned_uclause = uclause.substitute(subYX)
    friendsXX = Literal(Atom(friends, [X, X]), True)
    target_uclause = UnconstrainedClause([notlikesXa, friendsXX])
    assert(returned_uclause == target_uclause)

def test_cs():
    returned_cs = cs.substitute(subYX)
    target_cs = ConstraintSet([Xeqa, ~Xeqa, XinD])
    print(cs)
    print(returned_cs)
    assert(returned_cs == target_cs)
    assert(cs.is_satisfiable())
    assert(not returned_cs.is_satisfiable() )

def test_cclause():
    returned_cclause = cclause.substitute(subYX)
    friendsXX = Literal(Atom(friends, [X, X]), True)
    target_cclause = ConstrainedClause([notlikesXa, friendsXX], [X], cs.substitute(subYX))
    assert(returned_cclause == target_cclause)

def test_cnf():
    returned_cnf = cnf.substitute(subYX)
    target_cnf = CNF([uclause.substitute(subYX), cclause.substitute(subYX)])
    assert(returned_cnf == target_cnf)

