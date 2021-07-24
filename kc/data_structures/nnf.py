"""Classes for NNFs -- represented by their root nodes"""

from kc.data_structures import ConstrainedClause, UnconstrainedClause, Clause, Substitution, CNF

from abc import ABC, abstractmethod

from typing import Iterable, List, Set, Union, Dict, Tuple, Optional, FrozenSet
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import ConstraintSet, LogicalVariable, Literal, DomainVariable, ConstrainedAtom

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
    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'NNFNode'
        attributes['label'] = ''
        return attributes

    def _make_independent(self, c_atoms: List['ConstrainedAtom']) -> List['ConstrainedAtom']: 
        """Takes in a bunch of constrained atoms and returns an equivalent set that are all independent
        of each other.
        This is based on the Split algorithm, and is alluded to as a 'Modified Split' algorithm in the Smoothing
        section of the PhD"""
        viable_atom_pair = self._find_viable_atom_pair(c_atoms)
        if viable_atom_pair is None: # we are done if all are independent
            return c_atoms
        c_atom, other_c_atom = viable_atom_pair

        # if we have a viable atom pair, apply some preprocessing to the clauses to 
        # avoid variable name issues
        # DEBUG TODO: This is experimental
        c_atom = self._align_variables(c_atom, other_c_atom)
        cs1 = c_atom.cs
        cs2 = other_c_atom.cs

        mgu_eq_classes = c_atom.get_constrained_atom_mgu_eq_classes(other_c_atom)
        if mgu_eq_classes is None:
            raise ValueError(f"{c_atom = } and {other_c_atom  = } are independent but shouldn't be")

        theta = mgu_eq_classes.to_substitution()
        cs_theta = theta.to_constraint_set()
        joint_variables = c_atom.bound_vars.union(other_c_atom.bound_vars)

        return_clauses: List['ConstrainedAtom'] = []
        # loop over all constraints to negate
        for e in sorted(cs_theta.join(cs2)):
            # NOTE DEBUG: Trying - only include constraint if it is relevant to the clause
            if not e.terms[0] in c_atom.all_variables or e.terms[1] in c_atom.all_variables:
                continue
            not_e = ~e
            cs_rest = cs1.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            if cs_rest.is_satisfiable():
                return_clauses.append( ConstrainedAtom(c_atom.literals, joint_variables, cs_rest).propagate_equality_constraints() )
        return self._make_independent(return_clauses)

    def _find_viable_atom_pair(self, c_atoms: Iterable['ConstrainedAtom']
                              ) -> Optional[Tuple['ConstrainedAtom', 'ConstrainedAtom']]:
        for c_atom in c_atoms:
            for other_c_atom in c_atoms:
                if c_atom != other_c_atom and not c_atom.is_independent_from_other_clause(other_c_atom):
                    viable_atom_pair = (c_atom, other_c_atom)
                    return viable_atom_pair
        return None

    def _make_variables_different(self,
                                c_atom: 'ConstrainedAtom',
                                other_c_atom: 'ConstrainedAtom',
                                ) -> 'ConstrainedAtom':
        """MODIFIED FROM UNITPROP
        Apply preprocessing to the BOUND variables of two clauses to make them 
        all distinct. Also optionally apply this same substitution to another clause"""

        # all the variables that need to be substituted
        overlapping_variables: FrozenSet['LogicalVariable'] = c_atom.bound_vars.intersection(other_c_atom.bound_vars)

        old_other_c_atom = other_c_atom
        for variable in sorted(overlapping_variables):
            temp_cnf = CNF([c_atom, old_other_c_atom], names=None)  # taking advantage of existing methods in CNF
            sub_target = temp_cnf.get_new_logical_variable(variable.symbol)
            sub = Substitution([(variable, sub_target)])
            new_other_c_atom = old_other_c_atom.substitute(sub)
            if new_other_c_atom is None:
                raise ValueError('Substitution made unsatisfiable')
            else:
                old_other_c_atom = new_other_c_atom
        return old_other_c_atom

    def _align_variables(self,
                         c_atom: 'ConstrainedAtom',
                         other_c_atom: 'ConstrainedAtom'
                         ) -> 'ConstrainedAtom':
        """MODIFIED FROM UNITPROP
        'Line up' the variables in the  other_c_atom to match those in the first (c_atom), and
        apply this to the whole clause.
        This is done to make the mgu meaningful, rather than just having it rename variables to match.
        First, separate out all the bound variables.
        Then for each bound variable in the terms of this c_atom, we rename the other to match, as long as it 
        also is a bound variable."""
        other_c_atom = self._make_variables_different(c_atom, other_c_atom)
        old_other_c_atom = other_c_atom
        for term, other_term in zip(c_atom.atom.terms, other_c_atom.atom.terms):
            if term in c_atom.bound_vars and other_term in other_c_atom.bound_vars:
                assert(isinstance(other_term, LogicalVariable))  # hack for type checking
                sub = Substitution([(other_term, term)])
                new_other_c_atom = old_other_c_atom.substitute(sub)
                if new_other_c_atom is None:
                    raise ValueError("Shouldn't be unsatisfiable here")
                old_other_c_atom = new_other_c_atom
        return old_other_c_atom

    # def smooth(self, 

class TrueNode(NNFNode):
    """A class to represent True in the circuit. Needs no data"""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super(TrueNode, self).__init__([])

    def node_info(self) -> Dict[str, str]:
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

    def node_info(self) -> Dict[str, str]:
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

    def node_info(self) -> Dict[str, str]:
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

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'AndNode'
        and_string = '\u2227'
        attributes['label'] = and_string
        return attributes

class OrNode(ExtensionalNode):
    """A node representing an extensional OR operation."""

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'OrNode'
        or_string = '\u2228'
        attributes['label'] = or_string
        return attributes


class ForAllNode(IntensionalNode):
    """A node representing an intensional FORALL operation."""

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'ForAllNode'
        for_all_string = '\u2200'
        var_string = '{' + ', '.join(str(var) for var in self.bound_vars) + '}'
        attributes['label'] = f'{for_all_string}{var_string}, {self.cs}'
        return attributes


class ExistsNode(IntensionalNode):
    """A node representing an intensional EXISTS operation."""

    def node_info(self) -> Dict[str, str]:
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

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['type'] = 'EmptyNode'
        exists_string = '\u2203'
        attributes['label'] = ''
        return attributes
