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
import networkx as nx



# cython specifics -------------------------------------------------------------
import cython



# algorithms ###################################################################



# functions - encoding and decoding of graphs ##################################



# function: homogenizes nodes, node labels, and edgelabels of list of graphs ---
def encode_graphs(input_graphs = []):
    # description
    """
    > description: receives a list of inputs networkx (di)graphs, of the same type
    directed or undirected, and turns node names, node labels and edge labels into
    integers, considering their repetitions across the input list. All names and
    labels are encoded by integers starting from 0. The function returns a list with
    the copies of the input graphs preserving their order, together with dictionaries
    mapping integers into the original node_labels, edge_labels, and node_names. There
    are functions in networkx with similar behavior, but they are made for a given graph,
    while here we want to take into account all repetitions across the input list.

    > input:
    * input_graphs - list of networkx graphs or digraphs, all of the same type.

    > output:
    * encoded_graphs - list of transformed graphs.
    * node_name_encoding - dict from integers into original nodes (hashable objects).
    * node_label_encoding - dict from integers into original dicts of node labels.
    * edge_label_encoding - dict from integers into original dicts of edge labels.
    """
    # exception handling and input correctness
    test_list = [0, 0]
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    test_count_undir = 0
    test_count_dir = 0
    if(not type(input_graphs) in [type(test_list)]):
        raise(ValueError("gmapache: argument must be a list of networkx graphs or digraphs."))
    for test_entry in input_graphs:
        if(type(test_entry) not in [type(test_undir)]):
            if(type(test_entry) not in [type(test_dir)]):
                raise(ValueError("gmapache: elements in list must be networkx graphs or digraphs."))
            else:
                test_count_dir = test_count_dir + 1
        else:
            test_count_undir = test_count_undir + 1
    if(not test_count_undir == len(input_graphs)):
        if(not test_count_dir == len(input_graphs)):
            raise(ValueError("gmapache: elements in list must be networkx graphs or digraphs of the same type."))
    # output holders
    cdef list encoded_graphs = []
    cdef dict node_name_encoding = dict()    # from ints to node names
    cdef dict node_label_encoding = dict()   # from ints to node label-dicts
    cdef dict edge_label_encoding = dict()   # from ints to edge label-dicts
    # local variables (cython)
    cdef int i = 0
    cdef int N1 = 0
    cdef int N2 = 0
    cdef int encoded_node = 0
    cdef int encoded_node_a = 0
    cdef int encoded_node_b = 0
    cdef int encoded_label = 0
    cdef int new_node_label = 0
    cdef int new_edge_label = 0
    # local variables (python)
    cdef list all_nodes = []
    cdef list all_node_labels = []
    cdef list all_edge_labels = []
    cdef list node_label_encoding_inv = []   # from node label-dicts to ints (indices)
    cdef list edge_label_encoding_inv = []   # from edge label-dicts to ints (indices)
    cdef dict nodeInfo = dict()
    cdef dict edgeInfo = dict()
    cdef dict node_name_encoding_inv = dict()   # from node names to ints
    u = None
    v = None
    int_graph = None
    # homogenize node names across input graphs
    N1 = len(input_graphs)
    for i in range(0, N1):
        all_nodes = list(set(all_nodes + list(input_graphs[i].nodes())))
    N2 = len(all_nodes)
    for i in range(0, N2):
        node_name_encoding[i] = deepcopy(all_nodes[i])
        node_name_encoding_inv[deepcopy(all_nodes[i])] = i # nodes are hashable objects
    # homogenize node labels across input graphs
    for i in range(0, N1):
        for (v, nodeInfo) in list(input_graphs[i].nodes(data = True)):
            if(nodeInfo not in node_label_encoding_inv):
                node_label_encoding[new_node_label] = deepcopy(nodeInfo) # starting from 0
                node_label_encoding_inv.append(deepcopy(nodeInfo)) # labels are dictionaries
                new_node_label = new_node_label + 1
    # homogenize edge labels across input graphs
    for i in range(0, N1):
        for (u, v, edgeInfo) in list(input_graphs[i].edges(data = True)):
            if(edgeInfo not in edge_label_encoding_inv):
                edge_label_encoding[new_edge_label] = deepcopy(edgeInfo) # starting from 0
                edge_label_encoding_inv.append(deepcopy(edgeInfo)) # labels are dictionaries
                new_edge_label = new_edge_label + 1
    # create integerized copies of input graphs
    for i in range(0, N1):
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
                               GMNL = encoded_label) # gran_mapache_node_label
        # add encoded edges with encoded edge labels
        for (u, v, edgeInfo) in list(input_graphs[i].edges(data = True)):
            encoded_node_a = node_name_encoding_inv[u]
            encoded_node_b = node_name_encoding_inv[v]
            encoded_label = edge_label_encoding_inv.index(edgeInfo)
            int_graph.add_edge(encoded_node_a, encoded_node_b,
                               GMEL = encoded_label) # gran_mapache_edge_label
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
    > description: recieves encoded list of networkx (di)graphs, of the same type
    directed or undirected, and the dicts necessary for decodeding them. Returns
    the original graphs according to the decoding. This is basically the inverse
    operation of the function encode_graphs.

    > input:
    * encoded_graphs - list of transformed graphs.
    * node_name_encoding - dict from integers into original nodes (hashable objects).
    * node_label_encoding - dict from integers into original dicts of node labels.
    * edge_label_encoding - dict from integers into original dicts of edge labels.

    > output:
    * decoded_graphs - list of networkx graph or digraph objects.
    """
    # exception handling and input correctness
    test_list = [0, 0]
    test_dict = dict()
    test_count_undir = 0
    test_count_dir = 0
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    if(not type(encoded_graphs) in [type(test_list)]):
        raise(ValueError("gmapache: first argument must be a list of networkx graphs or digraphs encoded with granmapache."))
    if(not type(node_name_encoding) in [type(test_dict)]):
        raise(ValueError("gmapache: second argument must be a dictionary."))
    if(not type(node_label_encoding) in [type(test_dict)]):
        raise(ValueError("gmapache: third argument must be a dictionary."))
    if(not type(edge_label_encoding) in [type(test_dict)]):
        raise(ValueError("gmapache: fourth argument must be a dictionary."))
    for test_entry in encoded_graphs:
        if(type(test_entry) not in [type(test_undir)]):
            if(type(test_entry) not in [type(test_dir)]):
                raise(ValueError("gmapache: elements in list must be networkx graphs or digraphs."))
            else:
                test_count_dir = test_count_dir + 1
                for (test_node, test_info) in list(test_entry.nodes(data = True)):
                    if(test_node not in list(node_name_encoding.keys())):
                        raise(ValueError("gmapache: all the nodes of the input (di)graphs must be encoded by the input dictionaries."))
                    if("GMNL" not in list(test_info.keys())):
                        raise(ValueError("gmapache: one of the input (di)graphs is not encoded by granmapache."))
                    if(test_info["GMNL"] not in list(node_label_encoding.keys())):
                        raise(ValueError("gmapache: all the node-labels of the input (di)graphs must be encoded by the input dictionaries."))
                for (test_edge, test_info) in list(test_entry.edges(data = True)):
                    if("GMEL" not in list(test_info.keys())):
                        raise(ValueError("gmapache: one of the input (di)graphs is not encoded by granmapache."))
                    if(test_info["GMEL"] not in list(edge_label_encoding.keys())):
                        raise(ValueError("gmapache: all the edge-labels of the input (di)graphs must be encoded by the input dictionaries."))
        else:
            test_count_undir = test_count_undir + 1
            for (test_node, test_info) in list(test_entry.nodes(data = True)):
                if(test_node not in list(node_name_encoding.keys())):
                    raise(ValueError("gmapache: all the nodes of the input (di)graphs must be encoded by the input dictionaries."))
                if("GMNL" not in list(test_info.keys())):
                    raise(ValueError("gmapache: one of the input (di)graphs is not encoded by granmapache."))
                if(test_info["GMNL"] not in list(node_label_encoding.keys())):
                    raise(ValueError("gmapache: all the node-labels of the input (di)graphs must be encoded by the input dictionaries."))
            for (test_edge, test_info) in list(test_entry.edges(data = True)):
                if("GMEL" not in list(test_info.keys())):
                    raise(ValueError("gmapache: one of the input (di)graphs is not encoded by granmapache."))
                if(test_info["GMEL"] not in list(edge_label_encoding.keys())):
                    raise(ValueError("gmapache: all the edge-labels of the input (di)graphs must be encoded by the input dictionaries."))
    if(not test_count_undir == len(encoded_graphs)):
        if(not test_count_dir == len(encoded_graphs)):
            raise(ValueError("gmapache: elements in list must be networkx graphs or digraphs of the same type."))
    # output holders
    cdef list decoded_graphs = []
    # local variables (cython)
    cdef int i = 0
    cdef int u = 0
    cdef int v = 0
    cdef int N1 = 0
    cdef int decoded_node = 0
    cdef int decoded_node_a = 0
    cdef int decoded_node_b = 0
    cdef int decoded_label = 0
    # local variables (python)
    cdef dict nodeInfo = dict()
    cdef dict edgeInfo = dict()
    dec_graph = None
    # iterate decoding graphs
    N1 = len(encoded_graphs)
    for i in range(0, N1):
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
    # exception handling and input correctness
    test_list = [0, 0]
    test_tuple = (0, 0)
    test_dict = dict()
    if(not type(input_match) in [type(test_list)]):
        raise(ValueError("gmapache: first argument must be a list of 2-tuples."))
    if(not type(node_name_encoding) in [type(test_dict)]):
        raise(ValueError("gmapache: second argument must be a dictionary."))
    for test_entry in input_match:
        if(not type(test_entry) in [type(test_tuple)]):
            raise(ValueError("gmapache: all elements in input list must be tuples."))
        if(not len(test_entry) == 2):
            raise(ValueError("gmapache: all tuples in input list must be of lenght 2."))
        if(test_entry[0] not in list(node_name_encoding.values())):
            raise(ValueError("gmapache: input dictionary must encode all elements being matched."))
        if(test_entry[1] not in list(node_name_encoding.values())):
            raise(ValueError("gmapache: input dictionary must encode all elements being matched."))
    # output holders
    cdef list encoded_match = []
    # local variables (cython)
    cdef int i = 0
    cdef int N1 = 0
    cdef int index1 = 0
    cdef int index2 = 0
    # local variables (python)
    cdef dict node_name_encoding_inv = dict()
    # get encoding as array
    node_name_encoding_inv = {node_name_encoding[i]:i for i in range(len(node_name_encoding))}
    # encode match
    N1 = len(input_match)
    for i in range(0, N1):
        index1 = node_name_encoding_inv[input_match[i][0]]
        index2 = node_name_encoding_inv[input_match[i][1]]
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
    # exception handling and input correctness
    test_list = [0, 0]
    test_tuple = (0, 0)
    test_dict = dict()
    if(not type(encoded_match) in [type(test_list)]):
        raise(ValueError("gmapache: first argument must be a list of 2-tuples."))
    if(not type(node_name_encoding) in [type(test_dict)]):
        raise(ValueError("gmapache: second argument must be a dictionary."))
    for test_entry in encoded_match:
        if(not type(test_entry) in [type(test_tuple)]):
            raise(ValueError("gmapache: all elements in input list must be tuples."))
        if(not len(test_entry) == 2):
            raise(ValueError("gmapache: all tuples in input list must be of lenght 2."))
        if(test_entry[0] not in list(node_name_encoding.keys())):
            raise(ValueError("gmapache: input dictionary must encode all elements being matched."))
        if(test_entry[1] not in list(node_name_encoding.keys())):
            raise(ValueError("gmapache: input dictionary must encode all elements being matched."))
    # output holders
    cdef list decoded_match = []
    # local variables (cython)
    cdef int i = 0
    cdef int N1 = 0
    # local variables (python)
    node1 = None
    node2 = None
    # decode match
    N1 = len(encoded_match)
    for i in range(0, N1):
        node1 = node_name_encoding[encoded_match[i][0]]
        node2 = node_name_encoding[encoded_match[i][1]]
        decoded_match.append((deepcopy(node1), deepcopy(node2)))
    # end of function
    return(decoded_match)



################################################################################
################################################################################
