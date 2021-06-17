"""This file contains the abstract base class for compilation rules,
which the specific algorithms inherit from"""

from kc.data_structures import *
from abc import ABC, abstractmethod

class KCRule(ABC):
    """Abstract base class for compilation rules."""

    @abstractmethod
    def is_applicable(self, delta: 'CNF') -> bool:
        """Is this compilation rule applicable to the cnf 'delta'
        in its current state?"""
        pass

    @abstractmethod
    def apply(self, delta: 'CNF') -> 'NNFNode':
        """Apply this compilation rule to the cnf, returning an NNF"""
        pass

