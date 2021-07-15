"""This is a file for visualising the output of compilation.
I intend to use networkX with graphviz to produce NNF trees

Some inspiration is taken from https://stackoverflow.com/a/57512902/10005793"""

import matplotlib.pyplot as plt
import networkx as nx
import pydot
from networkx.drawing.nx_pydot import graphviz_layout

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kc.data_structures import NNFNode

from collections import deque

def build_nx_graph_from_nnf(root: 'NNFNode') -> 'DiGraph':
    """Take in an NNFNode and build a networkx graph for the circuit represented
    by the node. 
    The graph will be a tree, even though maybe it should be a more general DAG.
    The graph is constructed by a depth-first (I think, doesn't really matter) search of the children of the root"""

    graph = nx.DiGraph()
    unvisited_queue = deque([root])
    while len(unvisited_queue) > 0:
        current_node = unvisited_queue.pop()

        # add text for drawing
        node_data = current_node.node_info()
        node_string = f"{node_data['type']}\n{node_data['label']}"
        graph.add_node(current_node, label=node_string)

        for child in current_node.children:
            unvisited_queue.append(child)
            graph.add_edge(current_node, child)
    return graph
