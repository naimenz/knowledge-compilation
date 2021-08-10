"""Trying to compile a propositional example in two ways: 1-ary predicates with a constant, or 0-ary predicates with no constant"""
from kc.data_structures import *
from kc.compiler import Compiler
from kc.util import draw_nx_graph_from_nnf
from kc.parsing import write_nnf_to_txt

sun = Predicate('sun', 0)
rain = Predicate('rain', 0)
rainbow = Predicate('rainbow', 0)

Sun = Literal(Atom(sun, []))
Rain = Literal(Atom(rain, []))
Rainbow = Literal(Atom(rainbow, []))

clause = UnconstrainedClause([Rainbow, ~Sun, ~Rain])
cnf = CNF([clause])
compiler = Compiler()
nnf = compiler.compile(cnf)

# draw_nx_graph_from_nnf(nnf)
smoothed_nnf = nnf.do_smoothing(cnf)
# draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, "rainbow_theory")


