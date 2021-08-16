"""Attempting a WMI example to see if my SMT stuff works"""

from kc.data_structures import *
from kc.compiler import Compiler
from kc.util import draw_nx_graph_from_nnf

p = Predicate('p', 0)
x_15_20 = SMTPredicate('x', 0, 15, 20)
x_0_10 = SMTPredicate('x', 0, 0, 10)

p_literal = Literal(Atom(p, []))
x_15_20_literal = Literal(Atom(x_15_20, []))
x_0_10_literal = Literal(Atom(x_0_10, []))

clause = UnconstrainedClause([p_literal, x_0_10_literal])
# clause = UnconstrainedClause([x_15_20_literal, x_0_10_literal])

cnf = CNF([clause])

compiler = Compiler()

nnf = compiler.compile(cnf)
print(nnf)
# smoothed_nnf = nnf.get_smoothed_node()
smoothed_nnf = nnf.do_smoothing(cnf)

draw_nx_graph_from_nnf(nnf)
draw_nx_graph_from_nnf(smoothed_nnf)
