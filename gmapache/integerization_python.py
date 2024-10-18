################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: integerization                                                     #
#                                                                              #
# - Description: convert node "names", node attributes and edge attributes     #
#   into integers to simplify their comparison in the other rutines            #
#                                                                              #
################################################################################


# dependencies #################################################################


# already in python ------------------------------------------------------------
from copy import deepcopy


# not in python ----------------------------------------------------------------
import cython
import networkx as nx


# functions ####################################################################


# function: homogenizes nodes, node labels, and edgelabels of list of graphs ---
def encode_graphs(input_graphs = []):
    # description
    """
    - Receives a list of inputs graphs and turns each dictionary of
    node_labels and edge_labels into integers, considering their
    repetitions across the list.
    - Nodes are also renamed into integers starting from 1.
    - The function returns a list with the copies of the input_graph
    preserving their order, together with dictionaries from integers
    into the original node_labels, edge_labels, and node_names.
    """
    # output holders
    encoded_graphs = []
    node_name_encoding = dict()        # from ints to node names
    node_label_encoding = dict()       # from ints to node label-dicts
    edge_label_encoding = dict()       # from ints to edge label-dicts
    # cython variables
    cdef int i = 0
    # local variables
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


################################################################################
################################################################################
