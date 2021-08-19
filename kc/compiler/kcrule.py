"""This file contains the abstract base class for compilation rules,
which the specific algorithms inherit from"""

from kc.data_structures import *
from abc import ABC, abstractmethod

from typing import Tuple, Any, Optional
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.compiler import Compiler


class KCRule(ABC):
    """Abstract base class for compilation rules."""

    @classmethod
    @abstractmethod
    def is_applicable(cls, delta: 'CNF') -> Tuple[bool, Optional[Any]]:
        """Is this compilation rule applicable to the cnf 'delta'
        in its current state?
        Returns a boolean for whether it's applicable or not, and optional stored data"""
        pass

    @classmethod
    @abstractmethod
    def apply(cls, delta: 'CNF', stored_data: Any, compiler: 'Compiler') -> 'NNFNode':
        """Apply this compilation rule to the cnf, returning an NNF"""
        pass
