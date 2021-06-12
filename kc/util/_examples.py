"""Examples from the util files all in one place and not imported"""

from kc.data_structures import *
from kc.util import *

import random
random.seed(0)

# # utility.py
# print('utility.py')
# # preds = [Predicate('smokes', 2), Predicate('friends', 2), Predicate('fun', 1)]
# # constants: List['Constant'] = [Constant('a'), Constant('b'), Constant('c')]
# # variables = [LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')]

# # atoms = [Atom(preds[0], [variables[0], variables[1]]),
# #          Atom(preds[1], [variables[1], constants[0]]),
# #          Atom(preds[2], [constants[0]]),
# #          Atom(preds[1], [variables[0], variables[1]]),
# #          Atom(preds[1], [variables[1], variables[0]]),
# #         ]

# # literals = [Literal(atoms[0], True), Literal(atoms[1], False), Literal(atoms[2], False), Literal(atoms[3], True), Literal(atoms[4], True) ]
# # domains = [SetOfConstants(constants), SetOfConstants(constants[:2])]

# # constraints = [
# #                   # EqualityConstraint(variables[0], constants[0]),
# #                   InequalityConstraint(variables[1], constants[1]),
# #                   InequalityConstraint(variables[0], variables[1]),
# #                   InclusionConstraint(variables[0], domains[0]),
# #                   InclusionConstraint(variables[1], domains[0]),
# #                   # NotInclusionConstraint(variables[1], domains[1])
# #                   ]
# # constraint_set = ConstraintSet(constraints)
# # # print(constraint_set)
# # # print("get sol",get_solutions(constraint_set, variables[:2]))


# # constraints2 =[
# #                   # EqualityConstraint(variables[0], constants[0]),
# #                   InequalityConstraint(variables[1], constants[1]),
# #                   InequalityConstraint(variables[0], variables[1]),
# #                   InclusionConstraint(variables[0], domains[0]),
# #                   InclusionConstraint(variables[1], domains[0]),
# #                   # NotInclusionConstraint(variables[1], domains[1])
# #                   ]

# # constraint_set2 = ConstraintSet(constraints2)
# # # print(constraint_set2)
# # # print("get sol",get_solutions(constraint_set2, variables))


# # u_clause = UnconstrainedClause(literals[0:1] + literals[3:])
# # c_clause1 = ConstrainedClause(u_clause, variables[:2], constraint_set)
# # c_clause2 = ConstrainedClause(u_clause, variables[:2], constraint_set2)

# # c_atom1 = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], constraint_set)
# # c_atom2 = ConstrainedAtom(UnconstrainedClause([literals[-1]]), variables[:2], constraint_set2)
# # print("C_ATOMs:\n", c_atom1,'\n', c_atom2)
# # # print("ARGUMENTS:", get_solutions_to_constrained_atom(c_atom1))
# # print("Are the c-atoms independent?",constrained_atoms_independent(c_atom1, c_atom2))
# # print("Are the c-atoms subsumed?",constrained_atoms_subsumed(c_atom1, c_atom2))
# # print(get_constrained_atom_grounding(c_atom1))
# # print(get_constrained_atom_grounding(c_atom2))

# # # print(have_same_predicate(c_atoms[1], c_atoms[2]))
# # print("C_Clauses:\n", c_clause1,'\n',c_clause2)
# # print("Are the c-clauses independent?", constrained_clauses_independent(c_clause1, c_clause2))

# # subsumed = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], constraint_set)
# # subsumer = ConstrainedAtom(UnconstrainedClause([literals[-2]]), variables[:2], ConstraintSet(constraints[1:]))
# # print("Are the c-atoms subsumed?",constrained_atoms_subsumed(subsumer, subsumed))
# # print("Are the c-atoms subsumed?",constrained_atoms_subsumed(subsumed, subsumer))

# p_pred = Predicate('p', 3)
# a, b = Constant('a'), Constant('b')
# X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
# u_atom1 = Atom(p_pred, [X, X, a])
# u_atom2 = Atom(p_pred, [Y, Z, Z])

# mgu.py
print('mgu.py')
a, b, c = Constant('a'), Constant('b'), Constant('c')
X, Y, Z = LogicalVariable('X'), LogicalVariable('Y'), LogicalVariable('Z')
D, E = SetOfConstants([a, b, c]), SetOfConstants([b, c])

pred = Predicate('p', 3)
atom1 = Atom(pred, [X, X, a])
atom2 = Atom(pred, [Y, Y, a])

XeqY = EqualityConstraint(X, Y)
YeqX = EqualityConstraint(X, Y)
Yneqa = InequalityConstraint(Y, a)
Xeqa = EqualityConstraint(X, a)
Xeqb = EqualityConstraint(X, b)
Xeqc = EqualityConstraint(X, c)
XinD = InclusionConstraint(X, D)
YinD = InclusionConstraint(Y, D)
ZinD = InclusionConstraint(Z, D)
ZninE = NotInclusionConstraint(Z, E)
XinE = InclusionConstraint(X, E)
YinE = InclusionConstraint(Y, E)
ZinE = InclusionConstraint(Z, E)
XninE = NotInclusionConstraint(X, E)
cs1 = ConstraintSet([Yneqa, XinD, YinD, ZinD])
cs2 = ConstraintSet([~XeqY, Yneqa, XinD, YinD, ZinD])
# cs1 = ConstraintSet([XinD, ~Xeqa, ~Xeqb, ])
# cs2 = ConstraintSet([XinD, ])
c_atom1 = ConstrainedAtom(UnconstrainedClause([Literal(atom1, True)]), [X], cs1)
c_atom2 = ConstrainedAtom(UnconstrainedClause([Literal(atom2, True)]), [Y], cs2)

# print(c_atom1)
# print(c_atom2)
# print("CONSTRAINED")
# sub = get_constrained_atom_mgu_substitution(c_atom1, c_atom2)
# print(sub)
# if not sub is None:
#     print(substitution_to_constraint_set(sub))

print(is_satisfiable(ConstraintSet([XeqY, ~XeqY, YeqX, ~YeqX, XinE, YinE])))






