from kc.data_structures import *
from kc.compiler import ShatteredCompilation, Compiler

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
smokesalice = Literal(Atom(smokes, [alice]))
friendsXY = Literal(Atom(friends, [X, Y]))

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)
XeqY = EqualityConstraint(X, Y)
Xeqalice = InclusionConstraint(X, SetOfConstants([alice]))
Yeqalice = InclusionConstraint(Y, SetOfConstants([alice]))
cs = ConstraintSet([XinPeople, YinPeople])
# print(cs)
# print(cs.is_satisfiable())

clause1 = ConstrainedClause([~smokesX, ~friendsXY, smokesY], [X, Y], cs)
clause2 = ConstrainedClause([smokesalice], [], ConstraintSet([]))
cnf = CNF([clause1, clause2])
compiler = Compiler()
shattered_clauses_list = [ShatteredCompilation.shatter_clause(clause, terms, domains) for clause in cnf.clauses]
empty_set: Set['ConstrainedClause'] = set() # hack for type checking
flattened_shattered_clauses = empty_set.union(*shattered_clauses_list)
for clause in flattened_shattered_clauses:
    print(clause.cs)
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
