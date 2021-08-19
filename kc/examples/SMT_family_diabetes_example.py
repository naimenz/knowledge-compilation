"""This is the diabetes example from WFOMI (WITH SMT PREDICATES)"""
from kc.data_structures import *
from kc.compiler import *
from kc.util import build_nx_graph_from_nnf, draw_nx_graph_from_nnf
from kc.parsing import write_nnf_to_txt, make_auxiliary_predicates_for_clauses

X = LogicalVariable('X')
Y = LogicalVariable('Y')

alice = Constant('Alice')
bob = Constant('Bob')
charlie = Constant('Charlie')

People = RootDomain([alice, bob, charlie], 'person')

diabetes = Predicate('diabetes', 1)
family = Predicate('family', 2)
BMI = SMTPredicate('BMI', 1, 35, float('inf'))
# BMI = Predicate('BMI', 1)

diabetesX = Literal(Atom(diabetes, [X]))
diabetesalice = Literal(Atom(diabetes, [alice]))
BMIX = Literal(Atom(BMI, [X]))
BMIY = Literal(Atom(BMI, [Y]))
familyXY = Literal(Atom(family, [X, Y]))

XinPeople = InclusionConstraint(X, People)
YinPeople = InclusionConstraint(Y, People)

cs = ConstraintSet([XinPeople, YinPeople])

clause = ConstrainedClause([~BMIX, ~familyXY, BMIY], [X, Y], cs)
# print(clause)
auxiliary_clauses = make_auxiliary_predicates_for_clauses([clause])
print(auxiliary_clauses)
all_clauses = auxiliary_clauses + [UnconstrainedClause([diabetesalice])]
cnf = CNF(auxiliary_clauses)

# cnf.shattered = True  # hack for now because they don't seem to shatter in the PhD example
import time
tic = time.perf_counter()
compiler = Compiler()
nnf = compiler.compile(cnf)
# draw_nx_graph_from_nnf(nnf)

smoothed_nnf = nnf.do_smoothing(cnf)
toc = time.perf_counter()
print(f'Compiling to sda-DNNF took {toc - tic:0.4f} seconds')
# smoothed_nnf = nnf.get_smoothed_node()
draw_nx_graph_from_nnf(smoothed_nnf)

# in-place transformation of infs into bounds
minmax_dict = {'BMI': (10, 45)}
smoothed_nnf.replace_infinities(minmax_dict)
draw_nx_graph_from_nnf(smoothed_nnf)

write_nnf_to_txt(smoothed_nnf, 'diabetes_theory')
# write_nnf_and_weights_to_txt(smoothed_nnf, weights, 'diabetes_theory', 'diabetes_weights')
