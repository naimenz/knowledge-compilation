from kc.data_structures import *
from kc.compiler import ShatteredCompilation

X = LogicalVariable('X')
a, b, c = Constant('a'), Constant('b'), Constant('c')
alice, bob = Constant('alice'), Constant('bob')
Y, Z = LogicalVariable('Y'), LogicalVariable('Z')

D = SetOfConstants([a,b,c])
People = SetOfConstants([alice, bob])

free_vars = set([Y, Z])
consts = set([a, b])
terms = (free_vars, consts)
domains: set['DomainTerm'] = set([People, D])

# constraints = ShatteredCompilation.shatter_var(X, terms, domains)
# for constraint in constraints:
#     print(f'- {constraint}')
#     print(f'{constraint.is_satisfiable()=}')

Xeqa = InclusionConstraint(X, SetOfConstants([a]))
Xeqb = InclusionConstraint(X, SetOfConstants([b]))
XinD = InclusionConstraint(X, D)
XeqY = EqualityConstraint(X, Y)
XeqZ = EqualityConstraint(X, Z)
cs = ConstraintSet([~Xeqa, ~Xeqb, XinD, ~XeqY, ~XeqZ])
print(cs.is_satisfiable())
