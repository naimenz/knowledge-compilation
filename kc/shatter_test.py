from kc.data_structures import *
from kc.compiler import ShatteredCompilation

from typing import Set

X = LogicalVariable('X')
Y = LogicalVariable('Y')
alice = Constant('alice')
bob = Constant('bob')

# People = DomainVariable('People')
People = SetOfConstants([alice, bob])

free_vars: Set['LogicalVariable'] = set()
consts = set([alice])
terms = (free_vars, consts)
domains: Set['DomainTerm'] = set([People])

constraints = ShatteredCompilation.shatter_var(X, terms, domains)
# for constraint in constraints:
#     print(f'- {constraint}')

smokes = Predicate('smokes', 1)
friends = Predicate('friends', 2)
smokesX = Literal(Atom(smokes, [X]))
smokesY = Literal(Atom(smokes, [Y]))
friendsXY = Literal(Atom(friends, [X, Y]))

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
XeqY = EqualityConstraint(X, Y)
Xeqalice = InclusionConstraint(X, SetOfConstants([alice]))
Yeqalice = InclusionConstraint(Y, SetOfConstants([alice]))
cs = ConstraintSet([XinPeople, YinPeople, Xeqalice, Yeqalice, ~XeqY])
print(cs)
print(cs.is_satisfiable())

# clause1 = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs)
# clauses = ShatteredCompilation.shatter_clause(clause1, terms, domains)
# for clause in clauses:
#     clause.literals = []
#     if clause.cs.is_satisfiable():
#         print(f'- {clause}')
# Xeqa = InclusionConstraint(X, SetOfConstants([a]))
# Xeqb = InclusionConstraint(X, SetOfConstants([b]))
# XinD = InclusionConstraint(X, D)
# XeqY = EqualityConstraint(X, Y)
# XeqZ = EqualityConstraint(X, Z)
# cs = ConstraintSet([~Xeqa, ~Xeqb, XinD, ~XeqY, ~XeqZ])
