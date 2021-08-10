"""Testing out independence functions for smoothing"""
from kc.data_structures import *

X, Y = LogicalVariable('X'), LogicalVariable('Y')
a, b = Constant('a'), Constant('b')

friends = Predicate('friends', 2)
likes = Predicate('likes', 2)

friendsXY = Literal(Atom(friends, [X, Y]))
likesXY = Literal(Atom(likes, [X, Y]))

person = RootDomain([a, b], 'person')
person_top = DomainVariable('person-top', person, excluded_constants=[b])
person_bot = person_top.complement

Xinperson = InclusionConstraint(X, person)
Yinperson = InclusionConstraint(Y, person)

Xinperson_top = InclusionConstraint(X, person_top)
Xinperson_bot = InclusionConstraint(X, person_bot)
Yinperson_top = InclusionConstraint(Y, person_top)

c_atom = ConstrainedAtom([friendsXY], [X, Y], ConstraintSet([Xinperson, Yinperson]))
other_c_atom = ConstrainedAtom([friendsXY], [X, Y], ConstraintSet([Xinperson_top, Yinperson_top]))

likes_c_atom = ConstrainedAtom([likesXY], [X, Y], ConstraintSet([Xinperson, Yinperson_top]))
other_likes_c_atom = ConstrainedAtom([likesXY], [X, Y], ConstraintSet([Xinperson_bot, Yinperson_top]))

nnf = TrueNode()
subtracted = other_c_atom.subtract_from_c_atoms([c_atom, c_atom])
A = [c_atom, likes_c_atom]
B = [other_c_atom, other_likes_c_atom]
print(nnf.A_without_B(A, B))
