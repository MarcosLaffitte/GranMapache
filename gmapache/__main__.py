################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Description: suite for the analysis of bijections, morphisms, alignments   #
#   and other maps between graphs.                                             #
#                                                                              #
# - Example of usage                                                           #
#                                                                              #
# - Observations and ToDos:                                                    #
#   * for a large number of graphs, integerization with cython is kinda faster #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
import time
from copy import deepcopy


# not in python ----------------------------------------------------------------
import cython
import networkx as nx
import matplotlib.pyplot as plt


# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs
# from .common_subgraphs import bla
# from .verbosity import bli


# functions ####################################################################


# function: homogenizes nodes, node labels, and edgelabels of list of graphs ---
def encode_graphs_python(input_graphs = []):
    # description
    """
    > description: receives a list of inputs graphs and turns each dictionary
    of node_labels and edge_labels into integers, considering their repetitions
    across the list. Nodes are also renamed into integers starting from 1. The
    function returns a list with the copies of the input_graph preserving their
    order, together with dictionaries from integers into the original node_labels,
    edge_labels, and node_names.

    > input:
    * input_graphs - list of networkx graph or digraph objects (default to empty list).

    > output:
    * encoded_graphs - list of transformed graphs.
    * node_name_encoding - dict from integers into original nodes (hashable objects).
    * node_label_encoding - dict from integers into original dicts of node labels.
    * edge_label_encoding - dict from integers into original dicts of edge labels.
    """
    # output holders
    encoded_graphs = []
    node_name_encoding = dict()        # from ints to node names
    node_label_encoding = dict()       # from ints to node label-dicts
    edge_label_encoding = dict()       # from ints to edge label-dicts
    # local variables
    i = 0
    encoded_node = 0
    encoded_node_a = 0
    encoded_node_b = 0
    encoded_label = 0
    new_node_label = 1
    new_edge_label = 1
    all_nodes = []
    all_node_labels = []
    all_edge_labels = []
    node_name_encoding_inv = dict()    # from node names to ints
    node_label_encoding_inv = [None]   # from node label-dicts to ints (indices)
    edge_label_encoding_inv = [None]   # from edge label-dicts to ints (indices)
    int_graph = None
    # homogenize node names across input graphs
    for i in range(len(input_graphs)):
        all_nodes = list(set(all_nodes + list(input_graphs[i].nodes())))
    for i in range(len(all_nodes)):
        node_name_encoding[i+1] = deepcopy(all_nodes[i])
        node_name_encoding_inv[deepcopy(all_nodes[i])] = i+1
    # homogenize node labels across input graphs
    for i in range(len(input_graphs)):
        for (v, nodeInfo) in list(input_graphs[i].nodes(data = True)):
            if(nodeInfo not in node_label_encoding_inv):
                node_label_encoding[new_node_label] = deepcopy(nodeInfo)
                node_label_encoding_inv.append(deepcopy(nodeInfo))
                new_node_label = new_node_label + 1
    # homogenize edge labels across input graphs
    for i in range(len(input_graphs)):
        for (u, v, edgeInfo) in list(input_graphs[i].edges(data = True)):
            if(edgeInfo not in edge_label_encoding_inv):
                edge_label_encoding[new_edge_label] = deepcopy(edgeInfo)
                edge_label_encoding_inv.append(deepcopy(edgeInfo))
                new_edge_label = new_edge_label + 1
    # create integerized copies of input graphs
    for i in range(len(input_graphs)):
        # initialize graph with original type
        if(input_graphs[i].is_directed()):
            int_graph = nx.DiGraph()
        else:
            int_graph = nx.Graph()
        # add encoded nodes with encoded node labels
        for (v, nodeInfo) in list(input_graphs[i].nodes(data = True)):
            encoded_node = node_name_encoding_inv[v]
            encoded_label = node_label_encoding_inv.index(nodeInfo)
            int_graph.add_node(encoded_node,
                               GMNL = encoded_label)   # gran_mapache_node_label
        # add encoded edges with encoded edge labels
        for (u, v, edgeInfo) in list(input_graphs[i].edges(data = True)):
            encoded_node_a = node_name_encoding_inv[u]
            encoded_node_b = node_name_encoding_inv[v]
            encoded_label = edge_label_encoding_inv.index(edgeInfo)
            int_graph.add_edge(encoded_node_a, encoded_node_b,
                               GMEL = encoded_label)   # gran_mapache_edge_label
        # save graph in list
        encoded_graphs.append(deepcopy(int_graph))
    # end of function
    return(encoded_graphs, node_name_encoding, node_label_encoding, edge_label_encoding)


# function: decode graphs in list recovering original input graphs -------------
def decode_graphs_python(encoded_graphs = [],
                  node_name_encoding = dict(),
                  node_label_encoding = dict(),
                  edge_label_encoding = dict()):
    # description
    """
    > description: recieves encoded list of graphs and the dicts necessary for
    decodeding them. Returns the original graphs according to the decoding.
    This is basically the inverse operation of the function encode_graphs.

    > input:
    * encoded_graphs - list of transformed graphs.
    * node_name_encoding - dict from integers into original nodes (hashable objects).
    * node_label_encoding - dict from integers into original dicts of node labels.
    * edge_label_encoding - dict from integers into original dicts of edge labels.

    > output:
    * decoded_graphs - list of networkx graph or digraph objects.
    """
    # output holders
    decoded_graphs = []
    # local variables
    i = 0
    decoded_node = 0
    decoded_node_a = 0
    decoded_node_b = 0
    decoded_label = 0
    dec_graph = None
    # iterate decoding graphs
    for i in range(len(encoded_graphs)):
        # initialize graph with original type
        if(encoded_graphs[i].is_directed()):
            dec_graph = nx.DiGraph()
        else:
            dec_graph = nx.Graph()
        # compare vertices
        for (v, nodeInfo) in list(encoded_graphs[i].nodes(data = True)):
            decoded_node = deepcopy(node_name_encoding[v])
            decoded_label = deepcopy(node_label_encoding[nodeInfo["GMNL"]])   # gran_mapache_node_label
            dec_graph.add_nodes_from([(decoded_node, decoded_label)])
        # compare edges
        for (u, v, edgeInfo) in list(encoded_graphs[i].edges(data = True)):
            decoded_node_a = deepcopy(node_name_encoding[u])
            decoded_node_b = deepcopy(node_name_encoding[v])
            decoded_label = deepcopy(edge_label_encoding[edgeInfo["GMEL"]])   # gran_mapache_edge_label
            dec_graph.add_edges_from([(decoded_node_a, decoded_node_b, decoded_label)])
        # save graph in list
        decoded_graphs.append(deepcopy(dec_graph))
    # end of function
    return(decoded_graphs)


# main #########################################################################


if __name__ == "__main__":


    # buid first test graph
    G = nx.Graph()
    G.add_node("hola1", color = "blue", charge = -1, radius = 0.5)
    G.add_node("hola2", color = "red", charge = 1, radius = 0.75)
    G.add_node("hola3", color = "red", charge = 1, radius = 0.9)
    G.add_node("hola4", color = "blue", charge = -1, radius = 0.5)
    G.add_node("hola5", color = "green", charge = -1, radius = 0.25)
    G.add_edge("hola1", "hola2", weight = 10, bond = "single", strength = "weak")
    G.add_edge("hola2", "hola3", weight = 8, bond = "double", strength = "strong")
    G.add_edge("hola3", "hola4", weight = 1, bond = "triple", strength = "strong")
    G.add_edge("hola4", "hola5", weight = 10, bond = "single", strength = "weak")


    # buid second test graph
    H = nx.Graph()
    H.add_node("hola1", color = "red", charge = 1, radius = 0.5)
    H.add_node("hola6", color = "red", charge = 1, radius = 0.75)
    H.add_node("hola7", color = "red", charge = 1, radius = 0.9)
    H.add_node("hola8", color = "blue", charge = -1, radius = 0.5)
    H.add_node("hola9", color = "green", charge = -1, radius = 0.25)
    H.add_edge("hola1", "hola6", weight = 10, bond = "single", strength = "weak")
    H.add_edge("hola1", "hola7", weight = 7, bond = "triple", strength = "strong")
    H.add_edge("hola1", "hola8", weight = 10, bond = "triple", strength = "strong")
    H.add_edge("hola1", "hola9", weight = 3, bond = "single", strength = "strong")


    # build input graph
    input_list = [G, H]


    # testing encoding
    initial_time = time.time()
    graph_encoding = encode_graphs(input_list)
    final_time = time.time()
    time_encoding_cython = final_time - initial_time
    encoded_graph_list = graph_encoding[0]
    encoded_node_name = graph_encoding[1]
    encoded_node_label = graph_encoding[2]
    encoded_edge_label = graph_encoding[3]


    # testing decoding
    initial_time = time.time()
    recovered_graphs = decode_graphs(encoded_graph_list,
                                     encoded_node_name,
                                     encoded_node_label,
                                     encoded_edge_label)
    final_time = time.time()
    time_decoding_cython = final_time - initial_time


    # compare running times
    print("- Encoding Cython: ", time_encoding_cython)
    print("- Decoding Cython: ", time_decoding_cython)
    print("- Total Cython: ", time_encoding_cython + time_decoding_cython)


    # testing encoding
    initial_time = time.time()
    graph_encoding = encode_graphs_python(input_list)
    final_time = time.time()
    time_encoding_python = final_time - initial_time
    encoded_graph_list = graph_encoding[0]
    encoded_node_name = graph_encoding[1]
    encoded_node_label = graph_encoding[2]
    encoded_edge_label = graph_encoding[3]


    # testing decoding
    initial_time = time.time()
    recovered_graphs = decode_graphs_python(encoded_graph_list,
                                            encoded_node_name,
                                            encoded_node_label,
                                            encoded_edge_label)
    final_time = time.time()
    time_decoding_python = final_time - initial_time


    # compare running times
    print("- Encoding Python: ", time_encoding_python)
    print("- Decoding Python: ", time_decoding_python)
    print("- Total Python: ", time_encoding_python + time_decoding_python)


################################################################################
################################################################################
