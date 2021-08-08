"""Classes for NNFs -- represented by their root nodes"""

from kc.data_structures import ConstrainedClause, UnconstrainedClause, Clause, Substitution, CNF, LogicalVariable, ConstraintSet, ConstrainedAtom, Literal, LessThanConstraint, InequalityConstraint, SetOfConstants
from kc.util import get_element_of_set

from abc import ABC, abstractmethod

from typing import Iterable, List, Set, Union, Dict, Tuple, Optional, FrozenSet, TypeVar
from typing import TYPE_CHECKING

# to avoid circular imports that are just for type checking
if TYPE_CHECKING:
    from kc.data_structures import DomainVariable 

class NNFNode(ABC):
    """The abstract base class for all NNF nodes."""
    def __init__(self, children: Iterable['NNFNode'], from_smoothing=False) -> None:
        """TODO: work out the exact shape of this function"""
        self.children = tuple(children)
        self.parents: List['NNFNode'] = [] # only leaf nodes should have multiple parents
        # set parents of children?
        for child in self.children:
            child.parents.append(self)
        self.from_smoothing = from_smoothing  # a flag to keep track of whether this node is the result of smoothing

    @abstractmethod
    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = dict()
        attributes['smoothing'] = str(self.from_smoothing)
        return attributes

    @abstractmethod
    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        pass


    @abstractmethod
    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """Get the circuit atoms (as defined in my version of smoothing -- similar to atom_c in the PhD)
        for a given node, based on its children's circuit atoms"""
        pass


    @abstractmethod
    def get_smoothed_node(self) -> 'NNFNode':
        """Get a smoothed version of this node. 
        This recursively gets smoothed versions of this node's children."""
        pass

    def do_smoothing(self, cnf: 'CNF') -> 'NNFNode':
        """Initiate smoothing, which recursively calls get_smoothed_node.
        We pass in a cnf so we know which circuit atoms to aim for"""

        all_possible_atoms: Set['ConstrainedAtom'] = set.union(*(set(clause.get_constrained_atoms()) for clause in cnf.clauses))
        all_circuit_atoms = self.get_circuit_atoms()

        partially_smoothed_node = self.get_smoothed_node()
        # now add all atoms that were missed by the whole circuit
        missed_circuit_atoms = self.A_without_B(all_possible_atoms, all_circuit_atoms)
        smoothed_node = partially_smoothed_node.add_circuit_nodes(missed_circuit_atoms)
        return smoothed_node

    def add_circuit_nodes(self, circuit_nodes: Iterable['ConstrainedAtom']) -> 'NNFNode':
        """Add ForAll and Literal nodes for each ConstrainedAtom in the circuit_nodes.
        This is done to smooth OrNodes and ExistsNodes"""
        current_tail: 'NNFNode' = self
        for c_atom in sorted(circuit_nodes):
            true_literal_node = LiteralNode(get_element_of_set(c_atom.literals), from_smoothing=True)
            false_literal_node = LiteralNode(~get_element_of_set(c_atom.literals), from_smoothing=True)
            or_node = OrNode(true_literal_node, false_literal_node, from_smoothing=True)
            if len(c_atom.bound_vars) > 0 and c_atom.cs != ConstraintSet([]):
                smoothing_branch = ForAllNode(or_node, c_atom.bound_vars, c_atom.cs, from_smoothing=True)
                current_tail = AndNode(smoothing_branch, current_tail)
            else: 
                current_tail = AndNode(or_node, current_tail)
        return current_tail


    def _make_independent(self, c_atoms: Set['ConstrainedAtom']) -> Set['ConstrainedAtom']: 
        """Takes in a bunch of constrained atoms and returns an equivalent set that are all independent
        of each other.
        This is based on the Split algorithm, and is alluded to as a 'Modified Split' algorithm in the Smoothing
        section of the PhD"""
        # we are done if there is just one atom
        if len(c_atoms) == 1:
            return c_atoms
        viable_atom_pair = self._find_viable_atom_pair(c_atoms)
        if viable_atom_pair is None:  # we are done if all are independent 
            return c_atoms
        c_atom, other_c_atom = viable_atom_pair

        # if we have a viable atom pair, apply some preprocessing to the clauses to 
        # avoid variable name issues
        # DEBUG TODO: This is experimental
        other_c_atom = self._align_variables(c_atom, other_c_atom)
        cs = c_atom.cs
        other_cs = other_c_atom.cs

        mgu_eq_classes = c_atom.get_constrained_atom_mgu_eq_classes(other_c_atom)
        if mgu_eq_classes is None:
            raise ValueError(f"{c_atom = } and {other_c_atom  = } are independent but shouldn't be")

        theta = mgu_eq_classes.to_substitution()
        cs_theta = theta.to_constraint_set()
        joint_variables = c_atom.bound_vars.union(other_c_atom.bound_vars)

        return_clauses: Set['ConstrainedAtom'] = set()
        # loop over all constraints to negate
        for e in sorted(cs_theta.join(other_cs)):
            # NOTE DEBUG: Trying - only include constraint if it is relevant to the clause
            if not (e.terms[0] in c_atom.all_variables or e.terms[1] in c_atom.all_variables):
                continue
            not_e = ~e
            cs_rest = cs.join(ConstraintSet([not_e]))
            # before going further, check if the constraint set for this clause is even satisfiable
            if cs_rest.is_satisfiable():
                return_clauses.add( ConstrainedAtom(c_atom.literals, joint_variables, cs_rest).propagate_equality_constraints() )
        return self._make_independent(c_atoms.union(return_clauses) - set([c_atom]))

    def _find_viable_atom_pair(self, c_atoms: Iterable['ConstrainedAtom']
                              ) -> Optional[Tuple['ConstrainedAtom', 'ConstrainedAtom']]:
        """Find a pair of constrained atoms that are not independent"""
        for c_atom in sorted(c_atoms):
            for other_c_atom in sorted(c_atoms):
                if c_atom != other_c_atom and not c_atom.is_independent_from_other_clause(other_c_atom):
                    # find out if one subsumes another, so we know which atom to change
                    if c_atom.subsumes(other_c_atom):
                        viable_atom_pair = (c_atom, other_c_atom)
                    else:
                        viable_atom_pair = (other_c_atom, c_atom)
                    return viable_atom_pair
        return None

    # TODO: Rename this function
    def A_without_B(self, A: Set['ConstrainedAtom'], B: Set['ConstrainedAtom']) -> Set['ConstrainedAtom']:
        """Return a set of ConstrainedAtoms representing those in A not covered by B.
        The approach taken is to make everything from A independent or subsumed by B, then remove the subsumed parts"""
        split_A: Set['ConstrainedAtom'] = set()
        for c_atomA in sorted(A):
            independent_of_all = True
            if c_atomA in B:
                continue  # trivially covered if it is already in B
            else:
                for c_atomB in sorted(B):
                    if not c_atomA.is_independent_from_other_clause(c_atomB):
                        independent_of_all = False
                        if not c_atomB.subsumes(c_atomA):
                            # make c_atomA independent of c_atomB
                            independent_c_atomsA = self._make_independent(set([c_atomA, c_atomB])) - set([c_atomB])
                            split_A = split_A.union(independent_c_atomsA)
            if independent_of_all:
                split_A.add(c_atomA)
        # we do a second round of subsumes checks to avoid redundant newly created c_atoms
        # TODO: integrate this into first loop to make it more efficient
        new_A: Set['ConstrainedAtom'] = set()
        for c_atomA in sorted(split_A):
            subsumed = False
            for c_atomB in sorted(B):
                if c_atomB.subsumes(c_atomA):
                    subsumed = True
                    break
            if subsumed == False:
                new_A.add(c_atomA)
        # NOTE: the atoms aren't guaranteed to be independent when built this way, so we make them so
        return self._make_independent(new_A)


    def _find_ordered_viable_atom_pair(self, A: Iterable['ConstrainedAtom'], B: Iterable['ConstrainedAtom']
                              ) -> Optional[Tuple['ConstrainedAtom', 'ConstrainedAtom']]:
        """Find a pair of constrained atoms, with the first from A and the second from B,
        such that they are not not independent and the atom from B does not subsume the atom from A"""
        for c_atomA in sorted(A):
            for c_atomB in sorted(B):
                if c_atomA != c_atomB and not c_atomA.is_independent_from_other_clause(c_atomB):
                    if not c_atomB.subsumes(c_atomA):
                        viable_atom_pair = (c_atomA, c_atomB)
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


class TrueNode(NNFNode):
    """A class to represent True in the circuit. Needs no data"""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super(TrueNode, self).__init__([])

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """TrueNodes have no children and cover no constrained atoms."""
        return set()

    def get_smoothed_node(self) -> 'TrueNode':
        """TrueNodes do not change with smoothing"""
        return self

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(TrueNode, self).node_info()
        attributes['type'] = 'TrueNode'
        attributes['label'] = 'T'
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        raise NotImplementedError('WFOMI solver has no TrueNode')

class FalseNode(NNFNode):
    """A class to represent False in the circuit. Needs no data"""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super(FalseNode, self).__init__([])

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """FalseNodes have no children and cover no constrained atoms."""
        return set()

    def get_smoothed_node(self) -> 'FalseNode':
        """FalseNodes do not change with smoothing"""
        return self

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(FalseNode, self).node_info()
        attributes['type'] = 'FalseNode'
        attributes['label'] = 'F'
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        raise NotImplementedError('WFOMI solver has no FalseNode')

class LiteralNode(NNFNode):
    """A class to represent a single literal in the circuit. 
    Takes as input the literal it represents"""
    def __init__(self, literal: 'Literal', from_smoothing=False) -> None:
        super(LiteralNode, self).__init__([], from_smoothing=from_smoothing)
        self.literal = literal

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """LiteralNodes have no children and cover only their own literal.
        NOTE: we return it as a ConstrainedAtom despite it having no constraints
        NOTE IMPORTANT: We have to convert FreeVariables into non-free variables here
        so that they appear in bound_vars correctly when substituting"""
        c_atom = ConstrainedAtom([Literal(self.literal.atom, True)], [], ConstraintSet([]))
        # print(f'DEBUG: Literal node circuit_atoms:\n{c_atom}')
        return set([c_atom.replace_free_variables()])

    def get_smoothed_node(self) -> 'LiteralNode':
        """LiteralNodes do not change with smoothing"""
        return self

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(LiteralNode, self).node_info()
        attributes['type'] = 'LiteralNode'
        attributes['label'] = str(self.literal)
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        atom_string = str(self.literal.atom)
        # instead of ¬, WFOMI uses neg
        if self.literal.polarity == True:
            string = atom_string
        else:
            string = 'neg ' + atom_string
        return string


class ExtensionalNode(NNFNode):
    """Abstract subclass for AND and OR nodes
    In this implementation, ExtensionalNodes always have exactly two children,
    a left and a right."""
    def __init__(self, left: 'NNFNode', right: 'NNFNode', from_smoothing=False) -> None:
        super(ExtensionalNode, self).__init__((left, right), from_smoothing=from_smoothing)
        self.left = left
        self.right = right


class IntensionalNode(NNFNode):
    """Abstract subclass for FORALL and EXISTS nodes
    In this implementation, Intensional nodes always have exactly one child,
    and (possibly empty) bound variables and constraint sets"""
    child: 'NNFNode'

class AndNode(ExtensionalNode):
    """A node representing an extensional AND operation."""

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """Since Ands should be decomposable, the circuit atoms of the two branches are independent.
        Thus we can just get the circuit atoms of each branch and combine the sets"""
        left_circuit_atoms = self.left.get_circuit_atoms()
        right_circuit_atoms = self.right.get_circuit_atoms()
        all_circuit_atoms = left_circuit_atoms.union(right_circuit_atoms)
        # DEBUG TODO: checking that they really are independent
        # print(f"HERE:\n{self._make_independent(all_circuit_atoms) = }\n{                        all_circuit_atoms = }")
        all_circuit_atoms = self._make_independent(all_circuit_atoms)
        # assert(self._make_independent(all_circuit_atoms) == all_circuit_atoms)
        # print(f'DEBUG:AndNode circuit_atoms:\n{sorted(all_circuit_atoms)}')
        return all_circuit_atoms

    def get_smoothed_node(self) -> 'AndNode':
        """Smoothing AndNodes just makes a new AndNode with smoothed children"""
        return AndNode(self.left.get_smoothed_node(), self.right.get_smoothed_node())

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(AndNode, self).node_info()
        attributes['type'] = 'AndNode'
        and_string = '\u2227'
        attributes['label'] = and_string
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        string = 'and'  # all the hard work is done by the edges
        return string

class OrNode(ExtensionalNode):
    """A node representing an extensional OR operation."""

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """Since Ands should be decomposable, the circuit atoms of the two branches are independent.
        Thus we can just get the circuit atoms of each branch and combine the sets"""
        left_circuit_atoms = self.left.get_circuit_atoms()
        right_circuit_atoms = self.right.get_circuit_atoms()
        missing_from_left = self.A_without_B(right_circuit_atoms, left_circuit_atoms)
        missing_from_right = self.A_without_B(left_circuit_atoms, right_circuit_atoms)
        all_circuit_atoms = left_circuit_atoms.union(right_circuit_atoms)
        assert(self._make_independent(all_circuit_atoms) == all_circuit_atoms)
        # print(f'DEBUG:OrNode circuit_atoms:\n{all_circuit_atoms}')
        return all_circuit_atoms

    def get_smoothed_node(self) -> 'OrNode':
        """Smoothing OrNodes requires potentially adding AndNodes below each branch 
        with the missing circuit atoms"""
        left_circuit_atoms = self.left.get_circuit_atoms()
        right_circuit_atoms = self.right.get_circuit_atoms()
        missing_from_left = self.A_without_B(right_circuit_atoms, left_circuit_atoms)
        missing_from_right = self.A_without_B(left_circuit_atoms, right_circuit_atoms)
        smoothed_left = self.left.get_smoothed_node().add_circuit_nodes(missing_from_left)
        smoothed_right = self.right.get_smoothed_node().add_circuit_nodes(missing_from_right)
        return OrNode(smoothed_left, smoothed_right)

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(OrNode, self).node_info()
        attributes['type'] = 'OrNode'
        or_string = '\u2228'
        attributes['label'] = or_string
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        string = 'or'  # all the hard work is done by the edges
        return string

class ForAllNode(IntensionalNode):
    """A node representing an intensional FORALL operation."""
    def __init__(self, child: 'NNFNode', bound_vars: Iterable['LogicalVariable'], cs: 'ConstraintSet', from_smoothing: bool=False) -> None:
        super(ForAllNode, self).__init__((child,), from_smoothing=from_smoothing)
        self.child = child
        self.bound_vars = frozenset(bound_vars)
        self.cs = cs

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """For ForAll nodes, we just get the circuit atoms of the child and add on the constraint set and
        bound variables of this node"""
        child_circuit_atoms = self.child.get_circuit_atoms()
        # DEBUG TODO: This is a dirty hack where we replace '<' constraints with '!=' constraints because this should give us the same missing circuit atoms,
        # and the code is not set up to handle '<' constraints
        fixed_cs = ConstraintSet([c if not isinstance(c, LessThanConstraint) else InequalityConstraint(c.left_term, c.right_term) for c in self.cs.constraints])
        # due to how ISG/IPG build ForAllNodes, we have to replace the FreeVariables with LogicalVariables
        non_free_bound_vars = set(LogicalVariable(v.symbol) for v in self.bound_vars)
        circuit_atoms: Set['ConstrainedAtom'] = set(ConstrainedAtom(c.literals, c.bound_vars.union(non_free_bound_vars), c.cs.join(fixed_cs)).replace_free_variables() for c in child_circuit_atoms)
        # DEBUG TODO: checking that they really are independent
        assert(self._make_independent(circuit_atoms) == circuit_atoms)
        # print(f'DEBUG:ForAllNode circuit_atoms:\n{circuit_atoms}')
        return circuit_atoms

    def get_smoothed_node(self) -> 'ForAllNode':
        """Smoothing ForAllNodes just makes a new ForAllNode with smoothed child"""
        return ForAllNode(self.child.get_smoothed_node(), self.bound_vars, self.cs)

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(ForAllNode, self).node_info()
        attributes['type'] = 'ForAllNode'
        for_all_string = '\u2200'
        var_string = '{' + ', '.join(str(var) for var in self.bound_vars) + '}'
        attributes['label'] = f'{for_all_string}{var_string}, {self.cs}'
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        # TODO: Split double quantifiers into two for WFOMI
        if len(self.bound_vars) > 1:
            raise ValueError("WFOMI can't handle quantifying over two variables")
        else:
            bound_var = get_element_of_set(self.bound_vars)
            domain = self.cs.get_domain_for_variable(bound_var)
            # we need to include things that the variable is not equal to in the domain
            not_equal_term_strings = []
            for nc in self.cs.notinclusion_constraints:
                if nc.logical_term == bound_var and isinstance(nc.domain_term, SetOfConstants):
                    not_equal_term_strings.append(str(get_element_of_set(nc.domain_term.constants)))
            for ic in self.cs.inequality_constraints:
                if ic.left_term == bound_var:
                    not_equal_term_strings.append(str(ic.right_term))
                elif ic.right_term == bound_var:
                    not_equal_term_strings.append(str(ic.left_term))
            # we only need the slash and so on if there are any excluded terms
            if len(not_equal_term_strings) > 0:
                string = f"A{{{bound_var}}}{{{domain}/" + ",".join(not_equal_term_strings) + "}"
            else:
                string = f"A{{{bound_var}}}{{{domain}}}"
            return string

class ExistsNode(IntensionalNode):
    """A node representing an intensional EXISTS operation."""
    def __init__(self, child: 'NNFNode', bound_vars: Iterable['DomainVariable'], cs: 'ConstraintSet', from_smoothing: bool=False) -> None:
        super(ExistsNode, self).__init__((child,), from_smoothing=from_smoothing)
        self.child = child
        self.bound_vars = frozenset(bound_vars)
        self.cs = cs

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """For Exists nodes, we have to consider all the possible domains that could be quantified over. 
        TODO: For now I'm just substituting the parent domain back in for each c_atom, but I don't know if that's correct.
        """
        child_circuit_atoms = self.child.get_circuit_atoms()
        circuit_atoms: Set['ConstrainedAtom'] = self.substitute_parent_domain(child_circuit_atoms)
        # print(f'DEBUG:ExistsNode circuit_atoms:\n{circuit_atoms}')
        return circuit_atoms

    def get_smoothed_node(self) -> 'ExistsNode':
        """Smoothing AndNodes just makes a new AndNode with smoothed children"""
        # it's more efficient to redo this here than call self.get_circuit_atoms
        child_circuit_atoms = self.child.get_circuit_atoms()
        all_circuit_atoms: Set['ConstrainedAtom'] = self.substitute_parent_domain(child_circuit_atoms)

        # TODO: For now we loop to convergence, but there must be a quicker way
        missed_circuit_atoms = self.A_without_B(all_circuit_atoms, child_circuit_atoms)
        new_missed_circuit_atoms = self.A_without_B(missed_circuit_atoms, child_circuit_atoms)
        while missed_circuit_atoms != new_missed_circuit_atoms:
            missed_circuit_atoms = new_missed_circuit_atoms
            new_missed_circuit_atoms = self.A_without_B(missed_circuit_atoms, child_circuit_atoms)
        smoothed_child = self.child.get_smoothed_node().add_circuit_nodes(missed_circuit_atoms)

        return ExistsNode(smoothed_child, self.bound_vars, self.cs)

    def substitute_parent_domain(self, c_atoms: Iterable['ConstrainedAtom']) -> Set['ConstrainedAtom']:
        """Take a bunch of constrained atoms and, where applicable, substitute the quantified domain with its parent.
        NOTE TODO: For now, I also substitute the quantified domain's complement with the parent, but I'm not sure if this is appropriate"""
        subdomain = get_element_of_set(self.bound_vars)  # There should be a single bound domain variable
        parent_domain = get_element_of_set(self.cs.subset_constraints).right_term  # There should be only one here too
        new_c_atoms: Set['ConstrainedAtom'] = set()
        for c_atom in sorted(c_atoms):
            new_cs = c_atom.cs.substitute_domain_in_set_constraints(subdomain, parent_domain)
            new_cs = new_cs.substitute_domain_in_set_constraints(subdomain.complement, parent_domain)
            new_c_atom = ConstrainedAtom(c_atom.literals, c_atom.bound_vars, new_cs)
            new_c_atoms.add(new_c_atom)
        return new_c_atoms


    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(ExistsNode, self).node_info()
        attributes['type'] = 'ExistsNode'
        exists_string = '\u2203'
        var_string = '{' + ', '.join(str(var) for var in self.bound_vars) + '}'
        attributes['label'] = f'{exists_string}{var_string}, {self.cs}'
        return attributes

    def get_node_string(self) -> str:
        """Get the string that will be used to pass to the WFOMI computation"""
        # exists only ever have one bound var
        bound_var = get_element_of_set(self.bound_vars)
        parent_domain = get_element_of_set(self.cs.subset_constraints).right_term  # There should be only one here too

        string = f"E{{{bound_var}}}{{{parent_domain}}}"
        return string

class EmptyNode(NNFNode):
    """A node representing nothing. Realistically these shouldn't be 
    made, but I needed a placeholder.
    TODO: avoid these so I can delete the class again."""
    def __init__(self) -> None:
        # don't need any children but still want to have parents
        super().__init__([])

    def get_circuit_atoms(self) -> Set['ConstrainedAtom']:
        """EmptyNode is just a placeholder so this shouldn't even really be called."""
        return set()
        # raise NotImplementedError('Get rid of EmptyNodes')

    def get_smoothed_node(self) -> 'EmptyNode':
        """EmptyNode doesn't change"""
        return self
        # raise NotImplementedError('Get rid of EmptyNodes')

    def node_info(self) -> Dict[str, str]:
        """Get information about this node in a way that is used for drawing graphs"""
        attributes = super(EmptyNode, self).node_info()
        attributes['type'] = 'EmptyNode'
        exists_string = '\u2203'
        attributes['label'] = ''
        return attributes
