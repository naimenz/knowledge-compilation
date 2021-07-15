"""Classes for NNFs -- represented by their root nodes"""


from abc import ABC, abstractmethod

from typing import Iterable, List, Set, Union
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import ConstraintSet, LogicalVariable, Literal, DomainVariable

class NNFNode(ABC):
    """The abstract base class for all NNF nodes."""
    def __init__(self, children: Iterable['NNFNode']) -> None:
        """TODO: work out the exact shape of this function"""
        self.children = tuple(children)
        self.parents: List['NNFNode'] = [] # only leaf nodes should have multiple parents
        # set parents of children?
        for child in self.children:
            child.parents.append(self)

    @abstractmethod
    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'NNFNode'
        attributes['label'] = ''
        return attributes

class TrueNode(NNFNode):
    """A class to represent True in the circuit. Needs no data"""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super(TrueNode, self).__init__([])

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'TrueNode'
        attributes['label'] = 'T'
        return attributes

class FalseNode(NNFNode):
    """A class to represent False in the circuit. Needs no data"""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super(FalseNode, self).__init__([])

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'FalseNode'
        attributes['label'] = 'F'
        return attributes

class LiteralNode(NNFNode):
    """A class to represent a single literal in the circuit. 
    Takes as input the literal it represents"""
    def __init__(self, literal: 'Literal') -> None:
        super(LiteralNode, self).__init__([])
        self.literal = literal

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'LiteralNode'
        attributes['label'] = str(self.literal)
        return attributes


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
    def __init__(self, child: 'NNFNode', bound_vars: Iterable[Union['DomainVariable','LogicalVariable']], cs: 'ConstraintSet') -> None:
        super(IntensionalNode, self).__init__((child,))
        self.child = child
        self.bound_vars = frozenset(bound_vars)
        self.cs = cs

class AndNode(ExtensionalNode):
    """A node representing an extensional AND operation."""

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'AndNode'
        and_string = '\u2227'
        attributes['label'] = and_string
        return attributes

class OrNode(ExtensionalNode):
    """A node representing an extensional OR operation."""

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'OrNode'
        or_string = '\u2228'
        attributes['label'] = or_string
        return attributes


class ForAllNode(IntensionalNode):
    """A node representing an intensional FORALL operation."""

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'ForAllNode'
        for_all_string = '\u2200'
        var_string = '{' + ', '.join(str(var) for var in self.bound_vars) + '}'
        attributes['label'] = f'{for_all_string}{var_string}, {self.cs}'
        return attributes


class ExistsNode(IntensionalNode):
    """A node representing an intensional EXISTS operation."""

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'ExistsNode'
        exists_string = '\u2203'
        var_string = '{' + ', '.join(str(var) for var in self.bound_vars) + '}'
        attributes['label'] = f'{exists_string}{var_string}, {self.cs}'
        return attributes


class EmptyNode(NNFNode):
    """A node representing nothing. Realistically these shouldn't be 
    made, but I needed a placeholder.
    TODO: avoid these so I can delete the class again."""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super().__init__([])

    def node_info(self) -> str:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'EmptyNode'
        exists_string = '\u2203'
        attributes['label'] = ''
        return attributes
