"""Classes for NNFs -- represented by their root nodes"""

from abc import ABC, abstractmethod

from typing import Iterable, List

class NNFNode(ABC):
    """The abstract base class for all NNF nodes."""
    def __init__(self, children: Iterable['NNFNode']):
        """TODO: work out the exact shape of this function"""
        self.children = list(children)
        self.parents: List['NNFNode'] = [] # only leaf nodes should have multiple parents
        # set parents of children?
        for child in self.children:
            child.parents.append(self)
