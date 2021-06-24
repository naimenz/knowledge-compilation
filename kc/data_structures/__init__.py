from kc.data_structures.logicalterms import LogicalTerm, LogicalVariable, Constant
from kc.data_structures.domainterms import DomainTerm, SetOfConstants, DomainVariable

from kc.data_structures.cs_subs_eq_classes import Constraint, LogicalConstraint, SetConstraint, EqualityConstraint, InequalityConstraint, InclusionConstraint, NotInclusionConstraint, ConstraintSet
from kc.data_structures.cs_subs_eq_classes import Substitution
from kc.data_structures.cs_subs_eq_classes import EquivalenceClass, VariableEquivalenceClass, EquivalenceClasses

from kc.data_structures.literals import Literal, Atom, GroundAtom, Predicate
from kc.data_structures.clauses import Clause, UnconstrainedClause, ConstrainedClause, UnitClause, ConstrainedAtom
from kc.data_structures.nnf import NNFNode, ExtensionalNode, IntensionalNode, AndNode, OrNode, ForAllNode, ExistsNode
from kc.data_structures.cnf import CNF
