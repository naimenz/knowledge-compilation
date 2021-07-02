"""File for leaf construction. This isn't formally a compilation rule from the PhD,
but it is described in the Compile algorithm there. For consistency, we make it into a 
full rule."""

from kc.data_structures import *
from kc.compiler import KCRule
from kc.util import get_element_of_set

from typing import Tuple, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler

class LeafConstruction(KCRule):
    
    @classmethod
    def is_applicable(cls, cnf: 'CNF') -> Tuple[bool, Optional[str]]:
        """LeafConstruction is applicable if the theory is a single True, False, or literal.
        Returns True or False depending on the flag, plus None (no stored data needed)
        NOTE: I'm not quite sure what it means for a clause to be simply True
        TODO: Figure out the logic here instead of just returning False"""
        if len(cnf.clauses) == 1:
            # if it's a constrained clause, it can't be a single literal but it could be a tautology/contradiction
            clause: 'Clause' = get_element_of_set(cnf.clauses)
            if clause.is_contradiction():
                return True, "Contradiction"
            elif clause.is_tautology():
                return True, "Tautology"

            if isinstance(clause, UnconstrainedClause):
                if len(clause.literals) == 1:
                    return True, "Literal"
        return False, None

    @classmethod
    def apply(cls, cnf: 'CNF', node_type: str, compiler: 'Compiler') -> 'NNFNode':
        """Apply LeafConstruction and return an NNFNode"""
        if node_type == "Tautology":
            return TrueNode()
        elif node_type == "Contradiction":
            return FalseNode()
        elif node_type == "Literal":
            # must be a single unconstrained clause with one literal
            return LiteralNode(get_element_of_set(get_element_of_set(cnf.u_clauses).literals))
        else:
            raise ValueError(f"node_type should be one of 'Tautology', 'Contradiction' or 'Literal, not {node_type}")
        


