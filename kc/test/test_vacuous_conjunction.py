from kc.data_structures import *
from kc.compiler import VacuousConjunction, Compiler

class DummyCompiler(Compiler):
    """A dummy compiler so I can test the output without
    having to have all rules work"""
    def compile(self, cnf: 'CNF') -> 'CNF':
        return EmptyNode()

compiler = DummyCompiler()

X = LogicalVariable('X')
a = Constant('a')
b = Constant('b')
p = Predicate('p', 2)
pXa = Literal(Atom(p, [X, a]))

D = SetOfConstants([a, b])
XinD = InclusionConstraint(X, D)

cs = ConstraintSet([XinD])

c_clause = ConstrainedClause([pXa], set(), cs)
cnf = CNF([c_clause])

def test_is_applicable():
    applicable, stored_data = VacuousConjunction.is_applicable(cnf)
    assert(applicable == True)
    assert(stored_data is None)

# def test_apply():
#     returned_nnf = VacuousConjunction.apply(cnf, None, compiler)
#     u_clause = UnconstrainedClause([pXa])
#     child_nnf = NNFNode([CNF(u_clause)])
#     target_nnf = ForAllNode(EmptyNode(), set(), cs)

#     assert(returned_nnf == target_nnf)
