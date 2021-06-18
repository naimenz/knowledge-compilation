"""Examples from the util files all in one place and not imported"""

from kc.data_structures import *
from kc.util import *

# # utility.py
# print('utility.py')

# preds = [Predicate('smokes', 2), Predicate('friends', 2), Predicate('fun', 1)]
# constants: List['Constant'] = [Constant('a'), Constant('b'), Constant('c')]
# variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]

# atoms = [Atom(preds[0], [variables[0], variables[1]]),
#          Atom(preds[1], [variables[1], constants[0]]),
#          Atom(preds[2], [constants[0]]),
#          Atom(preds[1], [variables[0], variables[1]]),
#          Atom(preds[1], [variables[1], variables[0]]),
#         ]

# literals = [Literal(atoms[0], True), Literal(atoms[1], False), Literal(atoms[2], False), Literal(atoms[3], True), Literal(atoms[4], True) ]
# domains = [SetOfConstants(constants), SetOfConstants(constants[:2])]

# constraints = [
#                   # EqualityConstraint(variables[0], constants[0]),
#                   InequalityConstraint(variables[1], constants[1]),
#                   InequalityConstraint(variables[0], variables[1]),
#                   InclusionConstraint(variables[0], domains[0]),
#                   InclusionConstraint(variables[1], domains[0]),
#                   # NotInclusionConstraint(variables[1], domains[1])
#                   ]
# constraint_set = ConstraintSet(constraints)
# # print(constraint_set)
# # print("get sol",get_solutions(constraint_set, variables[:2]))


# constraints2 =[
#                   # EqualityConstraint(variables[0], constants[0]),
#                   InequalityConstraint(variables[1], constants[1]),
#                   InequalityConstraint(variables[0], variables[1]),
#                   InclusionConstraint(variables[0], domains[0]),
#                   InclusionConstraint(variables[1], domains[0]),
#                   # NotInclusionConstraint(variables[1], domains[1])
#                   ]

# constraint_set2 = ConstraintSet(constraints2)
# # print(constraint_set2)
# # print("get sol",get_solutions(constraint_set2, variables))


# u_clause = UnconstrainedClause(literals[0:1] + literals[3:])
# c_clause1 = ConstrainedClause(u_clause, variables[:2], constraint_set)
# c_clause2 = ConstrainedClause(u_clause, variables[:2], constraint_set2)

# c_atom1 = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], constraint_set)
# c_atom2 = ConstrainedAtom(UnconstrainedClause([literals[-1]]), variables[:2], constraint_set2)
# print("C_ATOMs:\n", c_atom1,'\n', c_atom2)
# # print("ARGUMENTS:", get_solutions_to_constrained_atom(c_atom1))
# print("Are the c-atoms independent?",constrained_atoms_independent(c_atom1, c_atom2))
# print("Are the c-atoms subsumed?",constrained_atoms_subsumed(c_atom1, c_atom2))
# print(get_constrained_atom_grounding(c_atom1))
# print(get_constrained_atom_grounding(c_atom2))

# # print(have_same_predicate(c_atoms[1], c_atoms[2]))
# print("C_Clauses:\n", c_clause1,'\n',c_clause2)
# print("Are the c-clauses independent?", constrained_clauses_independent(c_clause1, c_clause2))

# subsumed = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], constraint_set)
# subsumer = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], ConstraintSet(constraints[1:]))
# print("Are the c-atoms subsumed?",constrained_atoms_subsumed(subsumer, subsumed))
# print("Are the c-atoms subsumed?",constrained_atoms_subsumed(subsumed, subsumer))

# p_pred = Predicate('p', 3)
# a, b = Constant('a'), Constant('b')
# X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
# u_atom1 = Atom(p_pred, [X, X, a])
# u_atom2 = Atom(p_pred, [Y, Z, Z])

# # mgu.py
# print('mgu.py')
# a, b, c = Constant('a'), Constant('b'), Constant('c')
# X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
# D, E = SetOfConstants([a, b, c]), SetOfConstants([b, c])

# pred1 = Predicate('p', 3)
# pred2 = Predicate('q', 2)
# atom1 = Atom(pred1, [X, Y, a])
# atom2 = Atom(pred1, [Z, Z, a])
# atom3 = Atom(pred2, [X, Y])

# XeqY = EqualityConstraint(X, Y)
# YeqX = EqualityConstraint(X, Y)

# Yeqa = EqualityConstraint(Y, a)
# Yeqb = EqualityConstraint(Y, b)
# Yeqc = EqualityConstraint(Y, c)

# Xeqa = EqualityConstraint(X, a)
# Xeqb = EqualityConstraint(X, b)
# Xeqc = EqualityConstraint(X, c)

# XinD = InclusionConstraint(X, D)
# YinD = InclusionConstraint(Y, D)
# ZinD = InclusionConstraint(Z, D)

# XinE = InclusionConstraint(X, E)
# YinE = InclusionConstraint(Y, E)
# ZinE = InclusionConstraint(Z, E)

# cs1 = ConstraintSet([XinD, YinD, ZinD, YinE])
# cs2 = ConstraintSet([~Yeqa, XinD, YinD, ZinD])
# # cs1 = ConstraintSet([XinD, ~Xeqa, ~Xeqb, ])
# # cs2 = ConstraintSet([XinD, ])
# c_atom1 = ConstrainedAtom(UnconstrainedClause([Literal(atom1, True)]), [X], cs1)
# c_atom2 = ConstrainedAtom(UnconstrainedClause([Literal(atom2, True)]), [Z], cs2)

# print(c_atom1)
# print(c_atom2)
# print("CONSTRAINED")
# sub = get_constrained_atom_mgu_substitution(c_atom1, c_atom2)
# if not sub is None:
#     print("SATISFIABLE!")
#     print(substitution_to_constraint_set(sub))
# else:
#     print("NOT SATISFIABLE!")

# # splitting.py 
# print('splitting.py')
# cs1 = ConstraintSet([~Xeqc, XinD, YinD])
# cs2 = ConstraintSet([~Xeqb, XinD, YinD])
# atom = ConstrainedAtom(UnconstrainedClause([Literal(atom1, True)]), [X, Y], cs1)
# gamma = ConstrainedClause(UnconstrainedClause([Literal(atom1, False), Literal(atom3, True)]), [X, Y], cs2)
# print(gamma)
# print(split(gamma, atom))

# example 4.1
kiwi = Constant('kiwi')
penguin = Constant('penguin')
dog = Constant('dog')
pigeon = Constant('pigeon')

Animal = SetOfConstants([kiwi, penguin, pigeon, dog])
Bird = SetOfConstants([kiwi, penguin, pigeon])

X = LogicalVariable('X')
Y = LogicalVariable('Y')

flies = Predicate('flies', 1)
haswings = Predicate('haswings', 1)

flies_literalX = Literal(Atom(flies, [X]), True)
not_haswings_literalX = Literal(Atom(haswings, [X]), False)
flies_literalY = Literal(Atom(flies, [Y]), True)
not_haswings_literalY = Literal(Atom(haswings, [Y]), False)

Xeqkiwi = EqualityConstraint(X, kiwi)
Xeqpenguin = EqualityConstraint(X, penguin)
XinAnimal = InclusionConstraint(X, Animal)
XinBird = InclusionConstraint(X, Bird)

Yeqpenguin = EqualityConstraint(Y, penguin)
YinBird = InclusionConstraint(Y, Bird)

cs_gamma = ConstraintSet([~Xeqkiwi, XinAnimal])
# cs_a = ConstraintSet([~Yeqpenguin, YinBird])
cs_a = ConstraintSet([~Xeqpenguin, XinBird])

gamma = ConstrainedClause(UnconstrainedClause([~flies_literalX, not_haswings_literalX]), [X], cs_gamma)
atom = ConstrainedAtom(UnconstrainedClause([flies_literalX]), [X], cs_a)
# print("HERE WE GO")
# split_clauses = split(gamma, atom)
# for clause in split_clauses:
#     print(clause)

# example 4.2 (DOESN'T WORK BECAUSE OF FREE VARIABLES)
a, b, c = Constant('a'), Constant('b'), Constant('c')
X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
p, q = Predicate('p', 1), Predicate('q', 1)
D = SetOfConstants([a, b])
E = SetOfConstants([b, c])
Universe = SetOfConstants([a, b, c])

XinD = InclusionConstraint(X, D)
XinE = InclusionConstraint(X, E)
YeqZ = EqualityConstraint(Y, Z)
Xeqa = EqualityConstraint(X, a)
YinE = InclusionConstraint(Y, E)
XeqY = EqualityConstraint(X, Y)
YinUniverse = InclusionConstraint(Y, Universe)
ZinUniverse = InclusionConstraint(Z, Universe)

p_literal = Literal(Atom(p, [X]), True)
q_literal = Literal(Atom(q, [X]), True)

# cs_gamma = ConstraintSet([XinD, ~YeqZ, YinUniverse, ZinUniverse])
# cs_a = ConstraintSet([XinD, ~YeqZ, YinUniverse, ZinUniverse])

# gamma = ConstrainedClause(UnconstrainedClause([p_literal, q_literal]), [X], cs_gamma)
# atom = ConstrainedAtom(UnconstrainedClause([p_literal]), [X], cs_a)
# print(split(gamma, atom, Universe))

# another example
# cs_gammap = ConstraintSet([XinE, YinE, XeqY])
# cs_ap = ConstraintSet([XinD])

# gammap = ConstrainedClause(UnconstrainedClause([p_literal, q_literal]), [X], cs_gammap)
# atomp = ConstrainedAtom(UnconstrainedClause([p_literal]), [X], cs_ap)
# print(split(gammap, atomp, Universe))
# print(constrained_clauses_subsumed(atom, gamma))

# conditioning.py
# print('conditioning.py')
# clause = split_clauses[0]
# print(clause)
# print(atom)
# conditioned_clause = condition(clause, atom)
# print(conditioned_clause)

# unitprop.py
print('unitprop.py')
X = LogicalVariable('X')
Y = LogicalVariable('Y')
friends = Predicate('friends', 2)
likes = Predicate('likes', 2)
dislikes = Predicate('dislikes', 2)

a, b, c = Constant('a'), Constant('b'), Constant('c')
People = SetOfConstants([a, b])
Universe = SetOfConstants([a, b, c])

friendsXX = Literal(Atom(friends, [X, X]), True)
friendsXY = Literal(Atom(friends, [X, Y]), True)
likesXY = Literal(Atom(likes, [X, Y]), True)
dislikesXY = Literal(Atom(dislikes, [X, Y]), True)

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

clause1 = ConstrainedClause(UnconstrainedClause([friendsXY, dislikesXY]), [X, Y], ConstraintSet([XinPeople, YinPeople]))
clause2 = ConstrainedClause(UnconstrainedClause([~friendsXY, likesXY]), [X, Y], ConstraintSet([XinPeople, YinPeople]))
u = UnitClause(UnconstrainedClause([friendsXX]), [X], ConstraintSet([XinPeople]))

delta = CNF([clause1, clause2, u])

# clause1_split = split(clause1, get_constrained_atoms(u)[0])
# for c1 in clause1_split:
#     print(c1)
deltap = unitprop(delta, u, Universe)
print(deltap)

# forclift_independence.py
print('forclift_independence.py')
X = LogicalVariable('X')
Y = LogicalVariable('Y')
Z = LogicalVariable('Z')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')

friends = Predicate('friends', 2)
dislikes = Predicate('dislikes', 2)
likes = Predicate('likes', 2)

friendsXY = Literal(Atom(friends, [X, Y]), True)
dislikesXY = Literal(Atom(dislikes, [X, Y]), True)
friendsZZ = Literal(Atom(friends, [Z, Z]), True)
likesZZ = Literal(Atom(likes, [Z, Z]), True)
friendsYY = Literal(Atom(friends, [Y, Y]), True)
likesYY = Literal(Atom(likes, [Y, Y]), True)

uclause1 = UnconstrainedClause([friendsXY, dislikesXY])
uclause2 = UnconstrainedClause([~friendsZZ, likesZZ])
uclause3 = UnconstrainedClause([~friendsYY, likesYY])


People = SetOfConstants([alice, bob, charlie])
XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
ZinPeople = InclusionConstraint(Z, People)
XeqY = EqualityConstraint(X, Y)
Yeqalice = InclusionConstraint(Y, SetOfConstants([alice]))
Yeqbob = InclusionConstraint(Y, SetOfConstants([bob]))
Zeqbob = InclusionConstraint(Z, SetOfConstants([bob]))

cs1 = ConstraintSet([XinPeople, YinPeople, ~XeqY])
cs2 = ConstraintSet([ZinPeople])
cs3 = ConstraintSet([YinPeople, Yeqalice])
cs4 = ConstraintSet([ZinPeople, Zeqbob])

clause1 = ConstrainedClause(uclause1, [X, Y], cs1)
clause2 = ConstrainedClause(uclause2, [Z], cs2)
clause3 = ConstrainedClause(uclause3, [Y], cs3)
clause4 = ConstrainedClause(uclause2, [Z], cs4)

def run(cnf):
    res = tryIndependentSubtheories(cnf)
    if res is None:
        print(f"There are no independent subtheories in {cnf}")
    else:
        print("Recursing")
        run(res[0])
