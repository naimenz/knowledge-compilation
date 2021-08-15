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
        node_string = node_data['label']
        graph.add_node(current_node, label=node_string, smoothing=node_data['smoothing'])

        for child in current_node.children:
            unvisited_queue.append(child)
            graph.add_edge(current_node, child)
    return graph

def draw_nx_graph_from_nnf(root: 'NNFNode') -> None:
    graph = build_nx_graph_from_nnf(root)
    print(graph)

    pos = graphviz_layout(graph, prog="dot")

    label_pos = {}
    y_off = 10  # offset on the y axis

    for k, v in pos.items():
        label_pos[k] = (v[0], v[1]+y_off)
    node_labels = nx.get_node_attributes(graph, 'label')
    node_colours = ['blue' if n[1] == 'True' else 'green' for n in graph.nodes(data='smoothing')]
    nx.draw(graph, pos, node_color=node_colours)
    nx.draw_networkx_labels(graph, label_pos, labels=node_labels)
    plt.show()

def build_nx_graph_from_txt(file_name: 'str') -> 'DiGraph':
    with open(file_name + ".txt", 'r') as f:
        lines = f.readlines()
    node_lines = []
    edge_lines = []
    for line in lines:
        stripped_line = line.strip()
        if stripped_line[-1] == ';':
            edge_lines.append(stripped_line[:-1].split())
        else:
            node_lines.append(stripped_line.split(' '))

    graph = nx.DiGraph()
    for node_line in node_lines:
        node_name = node_line[0]
        node_data = "\n".join(node_line[1:])
        graph.add_node(node_name, label=node_name + '\n' + node_data)
    for edge_line in edge_lines:
        graph.add_edge(edge_line[0], edge_line[2])
    return graph

def draw_nx_graph_from_txt(file_name: 'str') -> None:
    """Just for testing purposes, draw a graph from a WFOMI txt file"""
    graph = build_nx_graph_from_txt(file_name)
    print(graph)

    pos = graphviz_layout(graph, prog="dot")

    label_pos = {}
    y_off = 10  # offset on the y axis

    for k, v in pos.items():
        label_pos[k] = (v[0], v[1]+y_off)
    node_labels = nx.get_node_attributes(graph, 'label')
    nx.draw(graph, pos, node_color='green')
    nx.draw_networkx_labels(graph, label_pos, labels=node_labels)
    plt.show()

if __name__ == "__main__":
    # draw_nx_graph_from_txt('../../WFOMI/solver/test_input/smokers/theory')
    draw_nx_graph_from_txt('../../WFOMI/solver/test_input/pipeline-smokers/query')
