"""This is a file to convert the NNF circuits produced by KC into 
the type of input expected by the WFOMI calculation code.

There are two ways this can be done -- by converting to a text file that is then parsed
by WFOMI, or by creating the WFOMI circuits directly.
We start with the former as it is likely simpler and makes debugging easier by keeping the two parts separate."""

from kc.data_structures import *
from kc.util import build_nx_graph_from_nnf
from kc.util import draw_nx_graph_from_nnf
from copy import copy, deepcopy

from typing import Optional, Tuple

from collections import deque

def nnf_to_wfomi_string(root: 'NNFNode') -> str:
    """Produce a string that, when written to a txt file, can be used as input
    for WFOMI"""
    nodes_string = ""
    edges_string = ""
    # doing a depth-first search to produce the node and edge strings
    global_index = 0  # keep track of which indices have already been used
    # get rid of things WFOMI can't handle
    root = prepare_nnf_for_wfomi(root)
    unvisited_queue = deque([(root, global_index)])
    global_index += 1
    while len(unvisited_queue) > 0:
        current_node, parent_index = unvisited_queue.pop()
        current_node_string = current_node.get_node_string()
        nodes_string += f"n{parent_index} {current_node_string}\n"

        for child in current_node.children:
            child_index = global_index
            # add an edge string for this connection
            edge_string = f"n{parent_index} -> n{child_index};\n"
            edges_string += edge_string
            unvisited_queue.append((child, child_index))
            global_index += 1

    string = nodes_string + edges_string
    return string

def prepare_nnf_for_wfomi(root: 'NNFNode') -> 'NNFNode':
    """There are some changes we must make to the circuit to make it fit
    with the WFOMI solver -- no TrueNodes, and no quantifiers over two bound variables.
    The approach here is to do a depth-first search until we encounter a problem, and then
    modify a copy of the circuit in-place"""
    # root = copy(root)  # hopefully this will also copy all children of root
    unvisited_queue = deque([root])
    while len(unvisited_queue) > 0:
        current_node = unvisited_queue.pop()
        # handling the case where we quantify over two variables
        if isinstance(current_node, ForAllNode) and len(current_node.bound_vars) == 2:
            returned_node: Optional['NNFNode']
            returned_node, current_node = replace_double_for_all(current_node)
            if returned_node is not None:
                root = returned_node

        # TODO: Make this less atrocious
        for child in current_node.children:
            # handling the specific case where this is an And of a TrueNode with something else
            # TODO maybe need to handle Or with True as well
            if isinstance(child, TrueNode) and isinstance(current_node, AndNode):
                returned_node = remove_true_node(child, current_node)
                if returned_node is not None:
                    root = returned_node
            else:
                unvisited_queue.append(child)
    return root

def remove_true_node(true_node: 'TrueNode', and_node: 'AndNode') -> Optional['NNFNode']:
    """Simplify the nnf IN PLACE by removing a redundant true/and node pair"""
    # if there are no parents, just return the other branch with no children
    if len(and_node.parents) == 0:
        if and_node.left == true_node:
            and_node.right.parents = []
            return and_node.right
        elif and_node.right == true_node:
            and_node.left.parents = []
            return and_node.left
    # otherwise we have to figure out what kind of node the parent is, and replace this AndNode and TrueNode with the other branch of the And.
    else:
        parent = and_node.parents[0]
        if isinstance(parent, IntensionalNode):
            if and_node.left == true_node:
                parent.child = and_node.right
                parent.children = (and_node.right,)
                and_node.right.parents = [parent]
            elif and_node.right == true_node:
                parent.child = and_node.left
                parent.children = (and_node.left,)
                and_node.left.parents = [parent]

        elif isinstance(parent, ExtensionalNode):
            if and_node.left == true_node:
                if parent.left == and_node:
                    parent.left = and_node.right
                    parent.children = (parent.left, parent.right)
                    and_node.right.parents = [parent]
                elif parent.right == and_node:
                    parent.right = and_node.left
                    parent.children = (parent.left, parent.right)
                    and_node.right.parents = [parent]

            elif and_node.right == true_node:
                if parent.left == and_node:
                    parent.left = and_node.left
                    parent.children = (parent.left, parent.right)
                    and_node.left.parents = [parent]
                elif parent.right == and_node:
                    parent.right = and_node.left
                    parent.children = (parent.left, parent.right)
                    and_node.left.parents = [parent]
    return None

def replace_double_for_all(for_all_node: 'ForAllNode') -> Tuple[Optional['ForAllNode'], 'ForAllNode'] :
    """Replace IN PLACE a for-all node over two variables with two for-all nodes
    If we have replaced the root, return that (or None). Additioally, return the new current node to search from
    NOTE TODO: This may cause problems with IPG"""
    bound_var, second_bound_var = sorted(for_all_node.bound_vars)
    domain, second_domain = for_all_node.cs.get_domain_for_variable(bound_var), for_all_node.cs.get_domain_for_variable(second_bound_var) 
    constraints = ConstraintSet([c for c in for_all_node.cs.set_constraints if c.logical_term == bound_var])
    # DEBUG TRying out putting logical constraints on larger domain
    second_constraints = ConstraintSet([c for c in for_all_node.cs.set_constraints if c.logical_term == second_bound_var])
    # the bigger domain gets the logical constraints
    if domain.is_strict_superset_of(second_domain):
        constraints = constraints.join(ConstraintSet(for_all_node.cs.logical_constraints))
    else:
        second_constraints = second_constraints.join(ConstraintSet(for_all_node.cs.logical_constraints))

    # second_constraints = ConstraintSet([c for c in for_all_node.cs.constraints if not (isinstance(c, SetConstraint) and c.logical_term == bound_var)])
    # build in reverse order so we can have the right children
    second_for_all = ForAllNode(for_all_node.child, [second_bound_var], second_constraints)
    for_all_node.child.parents = [second_for_all]
    new_for_all = ForAllNode(second_for_all, [bound_var], constraints)
    
    # update the children of the for_all_node
    for child in for_all_node.children:
        child.parents = [second_for_all]
    # now update the parent with the new for alls if there are any
    if len(for_all_node.parents) == 0:
        return new_for_all, second_for_all
    else:
        parent = for_all_node.parents[0]
        if isinstance(parent, IntensionalNode):
            parent.child, parent.children = new_for_all, (new_for_all,)
        elif isinstance(parent, ExtensionalNode):
            if parent.left == for_all_node:
                parent.left = new_for_all
                parent.children = (parent.left, parent.right)
            if parent.right == for_all_node:
                parent.right = new_for_all
                parent.children = (parent.left, parent.right)
    return None, second_for_all

def write_string_to_txt(string: str, file_name: str) -> None:
    """Write a string to a txt file"""
    with open(file_name + ".txt", 'w') as f:
        f.write(string)

def write_nnf_to_txt(root: 'NNFNode', file_name: str) -> None:
    string = nnf_to_wfomi_string(root)
    print(string)
    write_string_to_txt(string, file_name)

if __name__ == "__main__":
    X = LogicalVariable('X')
    Y = LogicalVariable('Y')
    People = RootDomain([Constant('alice')], 'People')
    clause1 = Literal(Atom(Predicate("hi", 1), [X]))
    clause2 = Literal(Atom(Predicate("bye", 1), [X]))
    root = AndNode(LiteralNode(clause1), AndNode(LiteralNode(clause2), TrueNode()))
    processed_root = prepare_nnf_for_wfomi(root)
    draw_nx_graph_from_nnf(root)
    draw_nx_graph_from_nnf(processed_root)

#     root1 = ForAllNode(AndNode(LiteralNode(clause2), TrueNode()), [LogicalVariable('X')], ConstraintSet([]))
#     processed_root1 = prepare_nnf_for_wfomi(root1)
#     draw_nx_graph_from_nnf(root1)
#     draw_nx_graph_from_nnf(processed_root1)
#     root2 = ForAllNode(AndNode(LiteralNode(clause2), LiteralNode(clause1)), [X, Y], ConstraintSet([InclusionConstraint(X, People), InclusionConstraint(Y, People), LessThanConstraint(X, Y)]))
#     processed_root2 = prepare_nnf_for_wfomi(root2)
#     draw_nx_graph_from_nnf(root2)
#     draw_nx_graph_from_nnf(processed_root2)

