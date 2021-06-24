"""Classes for NNFs -- represented by their root nodes"""


from abc import ABC, abstractmethod

from typing import Iterable, List, Set
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import ConstraintSet, LogicalVariable

class NNFNode(ABC):
    """The abstract base class for all NNF nodes."""
    def __init__(self, children: Iterable['NNFNode']) -> None:
        """TODO: work out the exact shape of this function"""
        self.children = tuple(children)
        self.parents: List['NNFNode'] = [] # only leaf nodes should have multiple parents
        # set parents of children?
        for child in self.children:
            child.parents.append(self)

class ExtensionalNode(NNFNode):
    """Abstract subclass for AND and OR nodes
    In this implementation, ExtensionalNodes always have exactly two children,
    a left and a right."""
    def __init__(self, left: 'NNFNode', right: 'NNFNode') -> None:
        super(ExtensionalNode, self).__init__((left, right))
        self.left = left
        self.right = right


class IntensionalNode(NNFNode):
    """Abstract subclass for FORALL and EXISTS nodes
    In this implementation, Intensional nodes always have exactly one child,
    and (possibly empty) bound variables and constraint sets"""
    def __init__(self, child: 'NNFNode', bound_vars: Set['LogicalVariable'], cs: 'ConstraintSet') -> None:
        super(IntensionalNode, self).__init__((child,))
        self.child = child
        self.bound_vars = bound_vars
        self.cs = cs

class AndNode(ExtensionalNode):
    """A node representing an extensional AND operation."""

class OrNode(ExtensionalNode):
    """A node representing an extensional OR operation."""

class ForAllNode(IntensionalNode):
    """A node representing an intensional FORALL operation."""

class ExistsNode(IntensionalNode):
    """A node representing an intensional EXISTS operation."""
