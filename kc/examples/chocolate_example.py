"""This is the belgian-likes-chocolate example (Ex 3.3, Eq 3.1)
In non-CNF form, this is forall X, X \in People: belgian(X) => likes(X, chocolate) """
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf
from kc.parsing import write_nnf_to_txt

X = LogicalVariable('X')

alice = Constant('alice')
bob = Constant('bob')
charlie = Constant('charlie')
chocolate = Constant('chocolate')

People = RootDomain([alice, bob, charlie], 'person')

likes = Predicate('likes', 2)
belgian = Predicate('belgian', 1)

likesXchocolate = Literal(Atom(likes, [X, chocolate]))
belgianX = Literal(Atom(belgian, [X]))


XinPeople = InclusionConstraint(X, People)

cs = ConstraintSet([XinPeople])

clause = ConstrainedClause([~belgianX, likesXchocolate], [X], cs)

cnf = CNF([clause])
# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
compiler = Compiler()
nnf = compiler.compile(cnf)
draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
# smoothed_nnf = nnf.get_smoothed_node()
draw_nx_graph_from_nnf(smoothed_nnf)
write_nnf_to_txt(smoothed_nnf, 'chocolate_theory')
