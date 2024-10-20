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



# algorithms ###################################################################



# functions - encoding and decoding of graphs ##################################



# function: homogenizes nodes, node labels, and edgelabels of list of graphs ---
def encode_graphs(input_graphs = []):
    # description
    """
    > description: receives a list of inputs networkx (di)graphs, of the same type directed
    or undirected, and turns each  dictionary of node_labels and edge_labels into integers,
    considering their repetitions across the list. Nodes are also renamed into integers
    starting from 1. The function returns a list with the copies of the input_graph
    preserving  their order, together with dictionaries from integers into the original
    node_labels, edge_labels, and node_names.

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
    nodeInfo = dict()
    edgeInfo = dict()
    node_name_encoding_inv = dict()    # from node names to ints
    node_label_encoding_inv = [None]   # from node label-dicts to ints (indices)
    edge_label_encoding_inv = [None]   # from edge label-dicts to ints (indices)
    u = None
    v = None
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
def decode_graphs(encoded_graphs = [],
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
    # cython variables
    cdef int i = 0
    cdef int u = 0
    cdef int v = 0
    # local variables
    decoded_node = 0
    decoded_node_a = 0
    decoded_node_b = 0
    decoded_label = 0
    nodeInfo = dict()
    edgeInfo = dict()
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
            decoded_node = node_name_encoding[v]
            decoded_label = node_label_encoding[nodeInfo["GMNL"]]   # gran_mapache_node_label
            dec_graph.add_nodes_from([(deepcopy(decoded_node), deepcopy(decoded_label))])
        # compare edges
        for (u, v, edgeInfo) in list(encoded_graphs[i].edges(data = True)):
            decoded_node_a = node_name_encoding[u]
            decoded_node_b = node_name_encoding[v]
            decoded_label = edge_label_encoding[edgeInfo["GMEL"]]   # gran_mapache_edge_label
            dec_graph.add_edges_from([(deepcopy(decoded_node_a), deepcopy(decoded_node_b), deepcopy(decoded_label))])
        # save graph in list
        decoded_graphs.append(deepcopy(dec_graph))
    # end of function
    return(decoded_graphs)



# functions - encoding and decoding of matches #################################



# function: encode match given an encoding of node names -----------------------
def encode_match(input_match = [],
                 node_name_encoding = dict()):
    # description
    """
    > description: receives a list of 2-tuples of nodes and a dictionary encoding
    nodes into integers, and returns the corresponding encoded list of 2-tuples.

    > input:
    * input_match - list of 2-tuples (x, y) representing a match between graphs.
    * node_name_encoding - dict from integers into the original nodes.

    > output:
    * encoded_match - list of 2-tuples of integers (i, j) encoding the nodes.
    """
    # output holders
    encoded_match = []
    # cython variables
    cdef int i = 0
    cdef int index1 = 0
    cdef int index2 = 0
    # local variables
    encoding_as_array = [None]
    # get encoding as array
    for i in range(1, len(node_name_encoding)+1):
        encoding_as_array.append(deepcopy(node_name_encoding[i]))
    # encode match
    for i in range(len(input_match)):
        index1 = encoding_as_array.index(input_match[i][0])
        index2 = encoding_as_array.index(input_match[i][1])
        encoded_match.append((index1, index2))
    # end of function
    return(encoded_match)



# function: decode match given an encoding of node names -----------------------
def decode_match(encoded_match = [],
                 node_name_encoding = dict()):
    # description
    """
    > description: receives a list of 2-tuples of integers and a dictionary mapping
    integers to nodes, and returns the corresponding decoded list of 2-tuples of nodes.

    > input:
    * encoded_match - list of 2-tuples (i, j) of integers encoding pairs of nodes.
    * node_name_encoding - dict from integers into the original nodes.

    > output:
    * decoded_match - list of 2-tuples (x, y) of nodes representing the match.
    """
    # output holders
    decoded_match = []
    # cython variables
    cdef int i = 0
    # local variables
    node1 = None
    node2 = None
    # decode match
    for i in range(len(encoded_match)):
        node1 = node_name_encoding[encoded_match[i][0]]
        node2 = node_name_encoding[encoded_match[i][1]]
        decoded_match.append((deepcopy(node1), deepcopy(node2)))
    # end of function
    return(decoded_match)



################################################################################
################################################################################
