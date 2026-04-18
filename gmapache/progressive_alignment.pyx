################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: progressive_alignment                                              #
#                                                                              #
# - Description: pairwise progressive alignemnt of labeld (di)graphs given to  #
#   the program inside a list of (labeled) networkx graphs or digraphs.        #
#                                                                              #
#  --------------------------------------------------------------------------  #
#                                                                              #
# - LICENSE:                                                                   #
#                                                                              #
#   This file is part of the academic work published in                        #
#            TBA                                                               #
#   and it is released under                                                   #
#            MIT License Copyright (c) 2024 Marcos E. González Laffitte        #
#   See LICENSE file in                                                        #
#            https://github.com/MarcosLaffitte/GranMapache                     #
#   for full license details.                                                  #
#                                                                              #
################################################################################





# C++ string encoding ##########################################################
# cython: c_string_type=unicode, c_string_encoding=utf8





# dependencies #################################################################





# already in python ------------------------------------------------------------
import random
from math import sqrt
from copy import deepcopy
from itertools import combinations





# not in python ----------------------------------------------------------------
import networkx as nx
import gklearn.kernels.structuralspKernel as sspk





# cython specifics -------------------------------------------------------------
import cython





# custom dependencies ----------------------------------------------------------
from .partial_maps import search_maximum_common_anchored_subgraphs
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match





# algorithms ###################################################################





# progressive graph alignment - wrapper ########################################





# function: --------------------------------------------------------------------
def anchored_progressive_graph_alignment(input_graphs = [],           # input list of length N of networkx graphs of the same type
                                         anchor_classes = [],         # non-empty list of non-empty N-tuples establishing anchor-equivalence-classes
                                         node_labels = False,         # consider node labels when evaluating pairwise alignments
                                         edge_labels = False,         # consider edge labels when evaluating pairwise alignments
                                         reachability = True,         # all nodes should be reachable from at least one anchor node
                                         ambiguous_edges = True,      # consider ambiguous edges when building alignments
                                         test_correctness = True,     # test correctness of the input
                                         verbose = True):             # set to False to omit printing messages of progress

    # description
    """> description: receives a list of networkx graphs or digraphs of
    the same type, and a list of anchor classes between them. The
    anchor (equivalence) classes are tuples of the same length of the
    input list of graphs, and the i-th entry of each tuple is a vertex
    of the i-th graph in the list. Then the progressive alignment of
    the graphs is carried while preserving these anchor classes in
    each pairwise alignment by computing an anchored-MCS containing
    the anchor classes or their corresponding images, i.e., each
    anchor vertex of an intermediate alignment corresponds to only one
    anchor class. IMPORTANT: the input anchor classes MUST be pairwise
    isomorphisms between induced subgraphs of the input graphs. If the
    test_correctness option is set to True all these N-1 isomorphisms
    will be tested before starting the alignment, by checking the
    adjacency and label preservation. If these do not hold the program
    will stop and raise an error. Additionally, if node_labels or
    edge_labels are set to False, the isomorphism test will only check
    adjacency but not the corresponding labels, but when doing a
    pairwise alignment of two input graphs or alignments of alignments,
    call them G and H, the program will overwrite the labels of H with
    the labels of G when building the alignment A(G, H) based on an
    internal order given to G and H. NOTE: if you are unsure whether
    the anchor classes respresent isomorphisms and do not want to used
    the option test_correctness for some reason, the easiest work around
    is to remove the edges they induce from the input graphs before
    calling this function. Then all anchored-MCS graphs will simply
    extend the matches between those empty induced subgraphs.

    > input:
    * input_graphs - list of length N of networkx (di)graphs of the same type
    * anchor_classes - non-empty list of N-tuples establishing equivalence-classes
    * node_labels - boolean for considering node labels when evaluating pairwise alignments
    * edge_labels - boolean for considering edge labels when evaluating pairwise alignments
    * reachability - boolean for considering all nodes should be reachable from at least
    one of the nodes in the anchor for each pairwise alignment
    * ambiguous_edges - boolean for considering ambiguous edges when building the alignments
    * test_correctness - boolean for testing correctness of the input
    * verbose - controls printing messages about intermediate outputs

    > output:
    * alignment_dictionary - python dictionary with the following nine string keys:
        "dendrogram" - dendrogram tuple representation of guide_tree_graph
                       NOTE: this tuple can be converted into a string simply by casting it with str()
        "input_graphs" - list of input graphs (leaves of alignment are the
                         indices of the graphs in this list)
        "alignment_graph" - final alignment graph
        "guide_tree_graph" - graph representation of the guide tree where the indices of input_graphs
                             are the leaves and inner nodes are tuples depicting (sub) dendrograms
                             NOTE: these tuples can be converted into a string simply by casting them with str()
        "consensus_graphs" - dict mapping strings "k/N" to node-induced subgraphs of
                             the final alignment for all k = 1 ,..., N
        "similarity_matrix" - dict mapping pairs (i, j) with i <= j, to the graph
                              kernel-similarity between the input graphs Gi and Gj
        "alignment_columns" - dict mapping nodes of alignment_graph to lists of length N with
                              entries equal to a node of the input graphs and None whenever
                              there is no correspondence
        "anchor_classes_map" - dict mapping each input anchor class to vertices of the final
                               alignment graph
        "intermediate_alignments" - dict mapping nodes of guide_tree_graph to tuples of the foloowing form
                                    (int. alignment graph, sub_alignment_columns, anchor_map, ambiguous neighbors, aligned_leaves)
    """

    # output holders (python)
    alignment_graph = None                          # final alignment graph
    guide_tree_graph = None                         # graph representation of the guide tree where the indices of input_graphs are leaves
    cdef dict alignment_dictionary = dict()         # dict saving the final alignment and guide tree graphs, plus the following 6 output holders and the input_graphs
    cdef tuple dendrogram = tuple()                 # dendrogram tuple representation of guide_tree_graph
    cdef dict similarity_matrix = dict()            # dict mapping pairs (i, j) with i <= j, to the graph kernel-similarity between the input graphs Gi and Gj
    cdef dict consensus_graphs = dict()             # dict mapping strings "k/N" to node-induced subgraphs of the final alignment
    cdef dict alignment_columns = dict()            # dict mapping nodes of alignment_graph to nodes of the input graphs and None whenever there is no correspondence
    cdef dict anchor_classes_map = dict()           # dict mapping each input anchor class to vertices of the final alignment graph
    cdef dict intermediate_alignments = dict()      # dict mapping nodes of guide_tree_graph to (int. alignment, sub_alignment_columns, anchor_map, ambiguous neighbors, aligned_leaves)

    # local variables (cython)
    cdef int i = 0
    cdef int k = 0
    cdef int N = 0
    cdef int order = 0
    cdef int total_anchor_classes = 0

    # local variables (python)
    results_bool = False
    input_correctness = False
    cdef int test_int = 0
    cdef tuple test_tuple = tuple()
    cdef tuple each_cherry = tuple()
    cdef list nodes = []
    cdef list anchor = []
    cdef list cherries = []
    cdef list results_MCSs = []
    cdef list anchored_MCS = []
    cdef list leaf_neighbors = []
    cdef list aligned_leaves = []
    cdef list aligned_leaves_G1 = []
    cdef list aligned_leaves_G2 = []
    cdef dict columns_G1 = dict()
    cdef dict columns_G2 = dict()
    cdef dict anchor_map = dict()               # dictionary mapping input anchor classes to vertices of each intermediate alignment graph
    cdef dict anchor_map_G1 = dict()
    cdef dict anchor_map_G2 = dict()
    cdef dict meta_anchor_map = dict()
    cdef dict ambiguous_neighbors = dict()
    cdef dict sub_alignment_columns = dict()    # dictionary mapping vertices of each alignment to list of vertices from input graphs
    cdef dict ambiguous_neighbors_G1 = dict()
    cdef dict ambiguous_neighbors_G2 = dict()
    cdef dict current_alignment_nodes = dict()
    n1 = None
    n2 = None
    G1 = None
    G2 = None
    each_node = None
    distance_graph = None
    current_alignment = None
    trimmed_guide_tree = None

    # test input correctness
    if(test_correctness):
        input_correctness = progressive_alignment_input_correctness(input_graphs, anchor_classes, node_labels, edge_labels, reachability, ambiguous_edges, verbose)
        if(not input_correctness):
            return(alignment_dictionary)
    else:
        # dummy value true just in case it can be used in the future for something else
        input_correctness = True

    # get number of input graphs and anchor classes
    N = len(input_graphs)
    total_anchor_classes = len(anchor_classes)

    # create clustering graph
    similarity_matrix, distance_graph = progressive_alignment_evaluate_graph_kernels(input_graphs, N, node_labels, edge_labels, verbose)

    # create guide tree
    guide_tree_graph = progressive_alignment_get_guide_tree(distance_graph, N)
    trimmed_guide_tree = deepcopy(guide_tree_graph)

    # prepare input graphs as (intermediate alignment) leaves of the guide tree following the input order
    for i in range(N):

        # reinitialize variables for i-th graph
        nodes = []
        anchor_map = dict()
        sub_alignment_columns = dict()

        # get order and vertices of i-th graph
        order = input_graphs[i].order()
        nodes = deepcopy(list(input_graphs[i].nodes(data = False)))

        # prepare trivial alignment columns
        for k in range(order):
            sub_alignment_columns[nodes[k]] = [None]*N
            sub_alignment_columns[nodes[k]][i] = deepcopy(nodes[k])

        # prepare trivial anchor map
        for k in range(total_anchor_classes):
            anchor_map[anchor_classes[k]] = anchor_classes[k][i]

        # initialize ambiguous neighbors place holder
        ambiguous_neighbors = dict()

        # build list of aligned leaves
        aligned_leaves = [i]

        # save input graph as intermediate graph
        intermediate_alignments[i] = (deepcopy(input_graphs[i]), deepcopy(sub_alignment_columns), deepcopy(anchor_map), deepcopy(ambiguous_neighbors), deepcopy(aligned_leaves))

    # progressive graph alignment based on guide-tree pruning
    while(trimmed_guide_tree.order() > 1):

        # reinitialize variables
        cherries = []

        # get order and vertices of i-th graph
        order = trimmed_guide_tree.order()
        nodes = deepcopy(list(trimmed_guide_tree.nodes(data = False)))

        # get cherries, i.e., pairs of leaves with the same parent
        for k in range(order):
            leaf_neighbors = [each_node for each_node in trimmed_guide_tree.neighbors(nodes[k]) if(trimmed_guide_tree.degree(each_node) == 1)]
            if(len(leaf_neighbors) == 2):
                cherries.append((nodes[k], leaf_neighbors[0], leaf_neighbors[1]))

        # pairwise-alignment of leaves in each cherry
        for each_cherry in cherries:

            # reinitialize variables
            anchor = []
            anchor_map = dict()
            meta_anchor_map = dict()
            sub_alignment_columns = dict()

            # unpack the graphs to be aligned
            G1 = intermediate_alignments[each_cherry[1]][0]
            G2 = intermediate_alignments[each_cherry[2]][0]

            # unpack alignment columns
            columns_G1 = intermediate_alignments[each_cherry[1]][1]
            columns_G2 = intermediate_alignments[each_cherry[2]][1]

            # unpack anchor maps
            anchor_map_G1 = intermediate_alignments[each_cherry[1]][2]
            anchor_map_G2 = intermediate_alignments[each_cherry[2]][2]

            # unpack ambiguous neighbors
            ambiguous_neighbors_G1 = intermediate_alignments[each_cherry[1]][3]
            ambiguous_neighbors_G2 = intermediate_alignments[each_cherry[2]][3]

            # unpack the (disjoint) lists of aligned leaves and concatenate them
            aligned_leaves_G1 = intermediate_alignments[each_cherry[1]][4]
            aligned_leaves_G2 = intermediate_alignments[each_cherry[2]][4]
            aligned_leaves = sorted(deepcopy(aligned_leaves_G1) + deepcopy(aligned_leaves_G2))

            # get anchor for alignment
            for k in range(total_anchor_classes):
                n1 = anchor_map_G1[anchor_classes[k]]
                n2 = anchor_map_G2[anchor_classes[k]]
                meta_anchor_map[anchor_classes[k]] = (n1, n2)
                anchor.append((n1, n2))

            # run MCS routine
            results_MCSs, results_bool = search_maximum_common_anchored_subgraphs(nx_G = G1, nx_H = G2, input_anchor = anchor,
                                                                                  node_labels = node_labels, edge_labels = edge_labels,
                                                                                  all_extensions = False, reachability = reachability,
                                                                                  ambiguous_neighbors_G = ambiguous_neighbors_G1,
                                                                                  ambiguous_neighbors_H = ambiguous_neighbors_G2)

            # unpack MCS
            anchored_MCS = results_MCSs[0]

            # build alignment graph
            current_alignment, current_alignment_nodes, comp_mcs_G1, comp_mcs_G2 = progressive_alignment_build_alignment(G1, G2, anchored_MCS)

            # build anchor map
            for k in range(total_anchor_classes):
                anchor_map[anchor_classes[k]] = current_alignment_nodes["MCS"][meta_anchor_map[anchor_classes[k]]]

            # build intermediate alignment columns
            sub_alignment_columns = progressive_alignment_get_subcolumns(N, current_alignment_nodes,
                                                                         anchored_MCS,
                                                                         comp_mcs_G1, comp_mcs_G2,
                                                                         columns_G1, columns_G2,
                                                                         aligned_leaves_G1, aligned_leaves_G2)

            # build ambiguous neighbors within the alignment only if necessary
            if(ambiguous_edges):
                ambiguous_neighbors = deepcopy(progressive_alignment_get_ambiguous_neighbors(current_alignment, sub_alignment_columns, aligned_leaves))
            else:
                ambiguous_neighbors = dict()

            # save alignment graph
            intermediate_alignments[each_cherry[0]] = (deepcopy(current_alignment),
                                                       deepcopy(sub_alignment_columns),
                                                       deepcopy(anchor_map),
                                                       deepcopy(ambiguous_neighbors),
                                                       deepcopy(aligned_leaves))

            # trim guide tree
            trimmed_guide_tree.remove_node(each_cherry[1])
            trimmed_guide_tree.remove_node(each_cherry[2])

    # unpack values for simplified look-ups after alignment
    dendrogram = deepcopy(list(trimmed_guide_tree.nodes())[0])
    alignment_graph = deepcopy(intermediate_alignments[dendrogram][0])
    alignment_columns = deepcopy(intermediate_alignments[dendrogram][1])
    anchor_classes_map = deepcopy(intermediate_alignments[dendrogram][2])

    # obtain consensus graphs
    consensus_graphs = progressive_alignment_consensus_graphs(alignment_graph, alignment_columns, N)

    # save output
    alignment_dictionary["dendrogram"] = dendrogram
    alignment_dictionary["input_graphs"] = input_graphs
    alignment_dictionary["alignment_graph"] = alignment_graph
    alignment_dictionary["guide_tree_graph"] = guide_tree_graph
    alignment_dictionary["consensus_graphs"] = consensus_graphs
    alignment_dictionary["similarity_matrix"] = similarity_matrix
    alignment_dictionary["alignment_columns"] = alignment_columns
    alignment_dictionary["anchor_classes_map"] = anchor_classes_map
    alignment_dictionary["intermediate_alignments"] = intermediate_alignments

    # end of function
    return(alignment_dictionary)





# functions - input correctness ################################################





# function: --------------------------------------------------------------------
def progressive_alignment_input_correctness(input_graphs, anchor_classes, node_labels, edge_labels, reachability, ambiguous_edges, verbose):

    # local variables (cython)
    cdef int i = 0

    # local variables (python)
    cdef tuple test_tuple = tuple()
    test_bool = False
    good_isos = False
    graphs_dir = False
    graphs_undir = False
    test_dir = nx.DiGraph()
    test_undir = nx.Graph()
    test_entry = None

    # check list of input graphs not emnpty
    if(len(input_graphs) == 0):
        raise(ValueError("gmapache: list of input graphs must be non-empty."))

    # check list of input graphs not emnpty
    if(len(input_graphs) == 1):
        raise(ValueError("gmapache: list of input graphs must contain at least two graphs."))

    # check consistency of list networkx graphs or digraphs (of the same type)
    for test_entry in input_graphs:
        # reinitialize variables
        graphs_dir = False
        graphs_undir = False
        # check undirected
        if(type(test_entry) in [type(test_undir)]):
            graphs_undir = True
        # check dir
        if(type(test_entry) in [type(test_dir)]):
            graphs_dir = True
        # consistency of being graph or digraph
        if((not graphs_undir) and (not graphs_dir)):
            raise(ValueError("gmapache: list of input graphs contains a non-networkx (di)graph."))
        # consistency of same type
        if(graphs_undir and graphs_dir):
            raise(ValueError("gmapache: list of input graphs must contain only Graphs or only DiGraphs from networkx."))

    # check list of anchor classes
    if(len(anchor_classes) == 0):
        raise(ValueError("gmapache: list of anchor classes must be non-empty."))

    # check consistency of anchor classes
    for test_entry in anchor_classes:
        # check entry is a tuple
        if(type(test_entry) not in [type(test_tuple)]):
            raise(ValueError("gmapache: list of anchor classes must contain only tuples."))
        # check length of entry
        if(not len(test_entry) == len(input_graphs)):
            raise(ValueError("gmapache: each anchor class must be of the same length as the list of input graphs."))
        # check consistency of tuples with each graph
        for i in range(len(test_entry)):
            if(test_entry[i] not in list(input_graphs[i].nodes())):
                raise(ValueError("gmapache: the i-th element of each anchor class must be a node in the i-th input graph."))

    # check consistency of boolean argument node_labels
    if(type(node_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument node_labels must be boolean."))

    # check consistency of boolean argument edge_labels
    if(type(edge_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument edge_labels must be boolean."))

    # check consistency of boolean argument reachability
    if(type(reachability) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument reachability must be boolean."))

    # check consistency of boolean argument ambiguous_edges
    if(type(ambiguous_edges) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument ambiguous_edges must be boolean."))

    # check consistency of boolean argument verbose
    if(type(verbose) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument verbose must be boolean."))

    # check consistency of the anchor classes as isomorphisms between induced subgraphs of the input graphs
    good_isos = progressive_alignment_test_anchor_classes(input_graphs, anchor_classes, node_labels, edge_labels)

    # end of function
    return(True)





# function: --------------------------------------------------------------------
def progressive_alignment_test_anchor_classes(input_graphs, anchor_classes, node_labels, edge_labels):

    # local variables (cython)
    cdef int i = 0
    cdef int k = 0
    cdef int C = 0
    cdef int T = 0

    # local variables (python)
    cdef list anchor_nodes_G_temp = []
    cdef list anchor_nodes_G_next = []
    cdef list anchor_edges_G_temp = []
    cdef list anchor_edges_G_next = []
    cdef dict forward_map = dict()
    cdef dict inverse_map = dict()
    cdef dict labeled_nodes_G_temp = dict()
    cdef dict labeled_nodes_G_next = dict()
    cdef dict labeled_edges_G_temp = dict()
    cdef dict labeled_edges_G_next = dict()
    u = None
    v = None
    labels = None
    G_temp = None
    G_next = None

    # get total of anchor classes
    C = len(anchor_classes)

    # get total of isomorphisms (pairs of graphs) to traverse
    T = len(input_graphs) - 1

    # get initial induced subgraph
    G_temp = deepcopy(input_graphs[0])
    anchor_nodes_G_temp = [anchor_classes[k][0] for k in range(C)]
    G_temp.remove_nodes_from([v for v in G_temp.nodes() if(not v in anchor_nodes_G_temp)])
    if(not G_temp.order() == C):
        raise(ValueError("gmapache: there is an anchor class with repeated elements, these must be isomorphisms."))

    # get containers for labeled test
    anchor_edges_G_temp = deepcopy(list(G_temp.edges(data = True)))
    if(node_labels):
        labeled_nodes_G_temp = {v : labels for (v, labels) in G_temp.nodes(data = True)}
    if(edge_labels):
        labeled_edges_G_temp = {(u, v) : labels for (u, v, labels) in anchor_edges_G_temp}

    # traverse input list taking consecutive graphs
    for i in range(T):

        # reinitialize map containers
        forward_map = dict()
        inverse_map = dict()

        # get next induced subgraph
        G_next = deepcopy(input_graphs[i+1])
        anchor_nodes_G_next = [anchor_classes[k][i+1] for k in range(C)]
        G_next.remove_nodes_from([v for v in G_next.nodes() if(not v in anchor_nodes_G_next)])
        if(not G_next.order() == C):
            raise(ValueError("gmapache: there is an anchor class with repeated elements, these must be isomorphisms."))

        # get containers for labeled test
        anchor_edges_G_next = deepcopy(list(G_next.edges(data = True)))
        if(node_labels):
            labeled_nodes_G_next = {v : labels for (v, labels) in G_next.nodes(data = True)}
        if(edge_labels):
            labeled_edges_G_next = {(u, v) : labels for (u, v, labels) in anchor_edges_G_next}

        # get forward and inverse maps (nx nodes must be hashable objects)
        for k in range(C):

            # get vertex and its image
            x = anchor_nodes_G_temp[k]
            u = anchor_nodes_G_next[k]

            # test consistency of node labels if necessary
            if(node_labels):
                if(not labeled_nodes_G_temp[x] == labeled_nodes_G_next[u]):
                    raise(ValueError("gmapache: node label preservation was requested, but not respected by the input anchors."))

            # save maps
            forward_map[x] = deepcopy(u)
            inverse_map[u] = deepcopy(x)

        # quick size test
        if(not G_temp.size() == G_next.size()):
            raise(ValueError("gmapache: one of the input anchors is not an isomorphism."))

        # test consistency of edges and edge labels in the forward direction
        for k in range(len(anchor_edges_G_temp)):

            # get vertices
            x = anchor_edges_G_temp[k][0]
            y = anchor_edges_G_temp[k][1]

            # get image vertices
            u = forward_map[x]
            v = forward_map[y]

            # test edge
            if(not G_next.has_edge(u, v)):
                raise(ValueError("gmapache: one of the input anchors does not preserve adjacency."))

            # test edge label if necessary
            if(edge_labels):
                if((u, v) in labeled_edges_G_next.keys()):
                    if(not labeled_edges_G_temp[(x, y)] == labeled_edges_G_next[(u, v)]):
                        raise(ValueError("gmapache: edge label preservation was requested, but not respected by the input anchors."))
                else:
                    if(not labeled_edges_G_temp[(x, y)] == labeled_edges_G_next[(v, u)]):
                        raise(ValueError("gmapache: edge label preservation was requested, but not respected by the input anchors."))

        # test consistency of edges and edge labels in the inverse direction
        for k in range(len(anchor_edges_G_next)):

            # get vertices
            u = anchor_edges_G_next[k][0]
            v = anchor_edges_G_next[k][1]

            # get preimage vertices
            x = inverse_map[u]
            y = inverse_map[v]

            # test edge
            if(not G_temp.has_edge(x, y)):
                raise(ValueError("gmapache: one of the input anchors does not preserve adjacency."))

            # test edge label if necessary
            if(edge_labels):
                if((x, y) in labeled_edges_G_temp.keys()):
                    if(not labeled_edges_G_temp[(x, y)] == labeled_edges_G_next[(u, v)]):
                        raise(ValueError("gmapache: edge label preservation was requested, but not respected by the input anchors."))
                else:
                    if(not labeled_edges_G_temp[(y, x)] == labeled_edges_G_next[(u, v)]):
                        raise(ValueError("gmapache: edge label preservation was requested, but not respected by the input anchors."))

        # update temp class
        G_temp = deepcopy(G_next)
        anchor_nodes_G_temp = deepcopy(anchor_nodes_G_next)
        anchor_edges_G_temp = deepcopy(anchor_edges_G_next)
        if(node_labels):
            labeled_nodes_G_temp = deepcopy(labeled_nodes_G_next)
        if(edge_labels):
            labeled_edges_G_temp = deepcopy(labeled_edges_G_next)

    # end of function
    return(True)





# functions - evaluation of graph kernels ######################################





# function: --------------------------------------------------------------------
def progressive_alignment_evaluate_graph_kernels(input_graphs, N, node_labels, edge_labels, verbose):

    # local variables (cython)
    cdef int k = 0

    # local variables (python)
    cdef list modified_nodes = []
    cdef list modified_edges = []
    cdef list modified_input_graphs = []
    cdef dict similarity_matrix = dict()
    self_score_Gi = 0
    self_score_Gj = 0
    distance_Gi_Gj = 0
    pair_score_Gi_Gj = 0
    sspk_matrix = None
    node_kernels = None
    edge_kernels = None
    modified_graph = None
    distance_graph = None

    # build appropiate encoded graphs
    if(node_labels and edge_labels):
        encoded_graphs, encoded_node_names, encoded_node_labels, encoded_edge_labels = encode_graphs(input_graphs)
    else:
        # traverse graphs and modified them accordingly
        for k in range(N):

            # turn off node labels if necessary
            modified_nodes = deepcopy(list(input_graphs[k].nodes(data = node_labels)))

            # turn off edge labels if necessary
            modified_edges = deepcopy(list(input_graphs[k].edges(data = edge_labels)))

            # create modified graph
            modified_graph = deepcopy(input_graphs[k])
            modified_graph.clear()
            modified_graph.add_nodes_from(deepcopy(modified_nodes))
            modified_graph.add_edges_from(deepcopy(modified_edges))

            # store modified graph
            modified_input_graphs.append(deepcopy(modified_graph))

        # encode modified graphs
        encoded_graphs, encoded_node_names, encoded_node_labels, encoded_edge_labels = encode_graphs(modified_input_graphs)

    # set kernel functions
    node_kernels = {"symb": kronecker_delta_one_vertices,
                    "nsymb": kronecker_delta_one_vertices,
                    "mix": kronecker_delta_two_vertices}
    edge_kernels = {"symb": kronecker_delta_one_edges,
                    "nsymb": kronecker_delta_one_edges,
                    "mix": kronecker_delta_two_edges}

    # evaluate kernel
    sspk_matrix = sspk.structuralspkernel(encoded_graphs,
                                          parallel = None,
                                          verbose = verbose,
                                          node_label = "GMNL",
                                          edge_label = "GMEL",
                                          node_kernels = node_kernels,
                                          edge_kernels = edge_kernels)[0]

    # build distance matrix and distance graph
    indices = list(range(N))
    distance_graph = nx.Graph()
    for (i, j) in combinations(indices, r = 2):

        # get distance
        self_score_Gi = sspk_matrix[i][i]
        self_score_Gj = sspk_matrix[j][j]
        pair_score_Gi_Gj = sspk_matrix[i][j]
        distance_Gi_Gj = sqrt(self_score_Gi + self_score_Gj - 2*pair_score_Gi_Gj)

        # build similarity matrix
        similarity_matrix[(i, j)] = sspk_matrix[i][j]

        # build distance from kernel-similarities
        distance_graph.add_edge(i, j, distance = distance_Gi_Gj)

    # end of function
    return(similarity_matrix, distance_graph)





# function: delta required by sspk for "symb" and "nsymb" for vertices ---------
def kronecker_delta_one_vertices(label1, label2):

    # end of function
    if(label1 == label2):
        return(1)
    else:
        return(0)





# function: delta required by sspk for "symb" and "nsymb" for edges ------------
def kronecker_delta_one_edges(label1, label2):

    # end of function
    if(label1 == label2):
        return(1)
    else:
        return(0)





# function: delta required by sspk for "mix" for vertices ----------------------
# NOTE: declared but not really used
def kronecker_delta_two_vertices(label1, label2, weight1, weight2):

    # end of function
    if((label1 == label2) and (weight1 == weight2)):
        return(1)
    else:
        return(0)





# function: delta required by sspk for "mix" for edges -------------------------
# NOTE: declared but not really used
def kronecker_delta_two_edges(label1, label2, weight1, weight2):

    # end of function
    if((label1 == label2) and (weight1 == weight2)):
        return(1)
    else:
        return(0)





# functions - build guide tree #################################################





# function: --------------------------------------------------------------------
def progressive_alignment_get_guide_tree(distance_graph, N):

    # local variables (cython)
    cdef int k = 0
    cdef int ds = 0
    cdef int mk = 0

    # local variables (python)
    test_int = 0
    new_distance = 0
    minimum_distance = 0
    cdef list other_nodes = []
    cdef list dendro_edges = []
    cdef tuple dendro_node = tuple()
    guide_tree = None
    dendro_graph = None

    # create guide tree and add leaves
    guide_tree = nx.Graph()
    guide_tree.add_nodes_from(list(range(N)))

    # initialize dendro graph to be modified
    dendro_graph = deepcopy(distance_graph)

    # iterative construction of dendrogram
    while(dendro_graph.order() > 1):

        # get size of dendrograph
        ds = dendro_graph.size()

        # get list of all edges
        dendro_edges = deepcopy(list(dendro_graph.edges(data = True)))

        # get arbitrary distance
        mk = 0
        minimum_distance = dendro_edges[0][2]["distance"]

        # traverse edges getting one with minimum distance
        for k in range(ds):

            # check if minimizing
            if(dendro_edges[k][2]["distance"] < minimum_distance):

                # update distance
                minimum_distance = dendro_edges[k][2]["distance"]

                # save index of current minimum
                mk = k

        # build dendro node
        if((type(dendro_edges[mk][0]) in [type(test_int)]) and (type(dendro_edges[mk][1]) in [type(test_int)])):
            if(dendro_edges[mk][0] < dendro_edges[mk][1]):
                dendro_node = (dendro_edges[mk][0], dendro_edges[mk][1])
            else:
                dendro_node = (dendro_edges[mk][1], dendro_edges[mk][0])
        else:
            dendro_node = (dendro_edges[mk][0], dendro_edges[mk][1])

        # get all the other nodes
        other_nodes = deepcopy(list(dendro_graph.nodes()))
        other_nodes.remove(dendro_node[0])
        other_nodes.remove(dendro_node[1])

        # contract dendro graph
        dendro_graph.add_node(dendro_node)
        for k in range(len(other_nodes)):
            new_distance = (dendro_graph[other_nodes[k]][dendro_node[0]]["distance"] + dendro_graph[other_nodes[k]][dendro_node[1]]["distance"])/2
            dendro_graph.add_edge(dendro_node, other_nodes[k], distance = new_distance)
        dendro_graph.remove_node(dendro_node[0])
        dendro_graph.remove_node(dendro_node[1])

        # add dendro_node to guide tree graph
        guide_tree.add_node(dendro_node)
        guide_tree.add_edge(dendro_node[0], dendro_node)
        guide_tree.add_edge(dendro_node[1], dendro_node)

    # end of function
    return(guide_tree)





# functions - build intermediate alignments ####################################





# function: --------------------------------------------------------------------
def progressive_alignment_build_alignment(G1, G2, anchored_MCS):

    # local variables (cython)
    cdef int k = 0
    cdef int A = 0
    cdef int N1 = 0
    cdef int N2 = 0
    cdef int M1 = 0
    cdef int M2 = 0

    # local variables (python)
    new_node_name = 0
    cdef tuple temp_pair = tuple()
    cdef tuple temp_node = tuple()
    cdef tuple temp_edge = tuple()
    cdef list nodes_G1 = []
    cdef list nodes_G2 = []
    cdef list edges_G1 = []
    cdef list edges_G2 = []
    cdef list side_mcs_G1 = []
    cdef list side_mcs_G2 = []
    cdef list comp_mcs_G1 = []
    cdef list comp_mcs_G2 = []
    cdef dict info = dict()
    cdef dict forward_map = dict()
    cdef dict inverse_map = dict()
    cdef dict current_alignment_nodes = dict()
    u = None
    v = None
    current_alignment = None

    # initialize current_alignment_nodes dicts
    current_alignment_nodes["MCS"] = dict()
    current_alignment_nodes["cG1"] = dict()
    current_alignment_nodes["cG2"] = dict()

    # get order of graphs and size of MCS
    A = len(anchored_MCS)
    N1 = G1.order()
    N2 = G2.order()
    M1 = G1.size()
    M2 = G2.size()

    # build base graph for alignment having the same
    # type as the input graphs, either G1 or G2
    current_alignment = deepcopy(G1)
    current_alignment.clear()

    # get nodes inside the MCS and asign a new name to them
    for k in range(A):
        # split sides
        side_mcs_G1.append(anchored_MCS[k][0])
        side_mcs_G2.append(anchored_MCS[k][1])
        # turn pairs into maps
        forward_map[anchored_MCS[k][0]] = anchored_MCS[k][1]
        inverse_map[anchored_MCS[k][1]] = anchored_MCS[k][0]
        # encode mcs nodes
        current_alignment_nodes["MCS"][anchored_MCS[k]] = new_node_name
        new_node_name = new_node_name + 1

    # get labeled nodes inside and outside the MCS wrt G1
    nodes_G1 = deepcopy(list(G1.nodes(data = True)))
    for k in range(N1):
        # test if node is or not in anchor
        if(nodes_G1[k][0] in side_mcs_G1):
            # get the (unique) pair maching the node
            temp_pair = (nodes_G1[k][0], forward_map[nodes_G1[k][0]])
            # get label node only if necessary
            if(not current_alignment.has_node(current_alignment_nodes["MCS"][temp_pair])):
                temp_node = (current_alignment_nodes["MCS"][temp_pair], nodes_G1[k][1])
                current_alignment.add_nodes_from([temp_node])
        else:
            # save complement of anchor in G1
            comp_mcs_G1.append(nodes_G1[k][0])
            # encode new node
            current_alignment_nodes["cG1"][nodes_G1[k][0]] = new_node_name
            temp_node = (new_node_name, nodes_G1[k][1])
            current_alignment.add_nodes_from([temp_node])
            new_node_name = new_node_name + 1

    # get labeled nodes inside and outside the MCS wrt G2
    nodes_G2 = deepcopy(list(G2.nodes(data = True)))
    for k in range(N2):
        # test if node is or not in anchor
        if(nodes_G2[k][0] in side_mcs_G2):
            # get the (unique) pair maching the node
            # NOTE: all pairs in the MCS should be already paired at this point, this extra
            # loop is only for sanity-check, and also when NOT preserving node labels the
            # ones from G2 will be overwritten by the ones in G1 for all nodes in the MCS
            temp_pair = (inverse_map[nodes_G2[k][0]], nodes_G2[k][0])
            # get label node only if necessary
            if(not current_alignment.has_node(current_alignment_nodes["MCS"][temp_pair])):
                temp_node = (current_alignment_nodes["MCS"][temp_pair], nodes_G2[k][1])
                current_alignment.add_nodes_from([temp_node])
        else:
            # save complement of anchor in G2
            comp_mcs_G2.append(nodes_G2[k][0])
            # encode new node
            current_alignment_nodes["cG2"][nodes_G2[k][0]] = new_node_name
            temp_node = (new_node_name, nodes_G2[k][1])
            current_alignment.add_nodes_from([temp_node])
            new_node_name = new_node_name + 1

    # get labeled edges inside and outside the MCS wrt G1
    edges_G1 = deepcopy(list(G1.edges(data = True)))
    for k in range(M1):
        if((edges_G1[k][0] in side_mcs_G1) and (edges_G1[k][1] in side_mcs_G1)):
            # map nodes
            u = current_alignment_nodes["MCS"][(edges_G1[k][0], forward_map[edges_G1[k][0]])]
            v = current_alignment_nodes["MCS"][(edges_G1[k][1], forward_map[edges_G1[k][1]])]
            # create edge if necessary, G1 overwrites G2
            if(not current_alignment.has_edge(u, v)):
                temp_edge = (u, v, edges_G1[k][2])
                current_alignment.add_edges_from([temp_edge])
        else:
            # map nodes if necessary
            if(edges_G1[k][0] in side_mcs_G1):
                u = current_alignment_nodes["MCS"][(edges_G1[k][0], forward_map[edges_G1[k][0]])]
                v = current_alignment_nodes["cG1"][edges_G1[k][1]]
                temp_edge = (u, v, edges_G1[k][2])
            else:
                if(edges_G1[k][1] in side_mcs_G1):
                    u = current_alignment_nodes["cG1"][edges_G1[k][0]]
                    v = current_alignment_nodes["MCS"][(edges_G1[k][1], forward_map[edges_G1[k][1]])]
                    temp_edge = (u, v, edges_G1[k][2])
                else:
                    u = current_alignment_nodes["cG1"][edges_G1[k][0]]
                    v = current_alignment_nodes["cG1"][edges_G1[k][1]]
                    temp_edge = (u, v, edges_G1[k][2])
            # create edge outside MCS
            if(not current_alignment.has_edge(u, v)):
                current_alignment.add_edges_from([temp_edge])

    # get labeled edges inside and outside the MCS wrt G2
    edges_G2 = deepcopy(list(G2.edges(data = True)))
    for k in range(M2):
        if((edges_G2[k][0] in side_mcs_G2) and (edges_G2[k][1] in side_mcs_G2)):
            # map nodes
            u = current_alignment_nodes["MCS"][(inverse_map[edges_G2[k][0]], edges_G2[k][0])]
            v = current_alignment_nodes["MCS"][(inverse_map[edges_G2[k][1]], edges_G2[k][1])]
            # create edge if necessary, G1 overwrites G2
            if(not current_alignment.has_edge(u, v)):
                temp_edge = (u, v, edges_G2[k][2])
                current_alignment.add_edges_from([temp_edge])
        else:
            # map nodes if necessary
            if(edges_G2[k][0] in side_mcs_G2):
                u = current_alignment_nodes["MCS"][(inverse_map[edges_G2[k][0]], edges_G2[k][0])]
                v = current_alignment_nodes["cG2"][edges_G2[k][1]]
                temp_edge = (u, v, edges_G2[k][2])
            else:
                if(edges_G2[k][1] in side_mcs_G2):
                    u = current_alignment_nodes["cG2"][edges_G2[k][0]]
                    v = current_alignment_nodes["MCS"][(inverse_map[edges_G2[k][1]], edges_G2[k][1])]
                    temp_edge = (u, v, edges_G2[k][2])
                else:
                    u = current_alignment_nodes["cG2"][edges_G2[k][0]]
                    v = current_alignment_nodes["cG2"][edges_G2[k][1]]
                    temp_edge = (u, v, edges_G2[k][2])
            # create edge outside MCS
            if(not current_alignment.has_edge(u, v)):
                current_alignment.add_edges_from([temp_edge])

    # end of function
    return(current_alignment, current_alignment_nodes, comp_mcs_G1, comp_mcs_G2)





# functions - sub alignment columns ############################################





# function: --------------------------------------------------------------------
def progressive_alignment_get_subcolumns(N, current_alignment_nodes, anchored_MCS, comp_mcs_G1, comp_mcs_G2, columns_G1, columns_G2, aligned_leaves_G1, aligned_leaves_G2):

    # local variables (cython)
    cdef int k = 0
    cdef int aM = 0
    cdef int c1 = 0
    cdef int c2 = 0
    cdef int l1 = 0
    cdef int l2 = 0
    cdef int entry = 0

    # local variables (python)
    cdef list temp_list = []
    cdef dict sub_alignment_columns = dict()     # dictionary mapping vertices of each alignment to list of vertices from input graphs
    n1 = None
    n2 = None

    # get cardinalities of mcs, its complements in each graph and aligned leaves
    aM = len(anchored_MCS)
    c1 = len(comp_mcs_G1)
    c2 = len(comp_mcs_G2)
    l1 = len(aligned_leaves_G1)
    l2 = len(aligned_leaves_G2)

    # build alignment columns for the vertices in the anchored MCS
    for k in range(aM):
        # initialize empty alignment column
        temp_list = [None]*N
        # split nodes in the pair
        n1 = anchored_MCS[k][0]
        n2 = anchored_MCS[k][1]
        # traverse valid entries in the alignment columns wrt G1 (disjoint with G2)
        for i in range(l1):
            # get entry of column being checked
            entry = aligned_leaves_G1[i]
            # inherit entry
            temp_list[entry] = columns_G1[n1][entry]
        # traverse valid entries in the alignment columns wrt G2 (disjoint with G1)
        for i in range(l2):
            # get entry of column being checked
            entry = aligned_leaves_G2[i]
            # inherit  entry
            temp_list[entry] = columns_G2[n2][entry]
        # save alignment column of the new node
        sub_alignment_columns[current_alignment_nodes["MCS"][(n1, n2)]] = deepcopy(temp_list)

    # build alignment columns for the vertices in the complement of the MCS in G1
    for k in range(c1):
        # get node from G1
        n1 = comp_mcs_G1[k]
        # copy alignment column
        sub_alignment_columns[current_alignment_nodes["cG1"][n1]] = deepcopy(columns_G1[n1])

    # build alignment columns for the vertices in the complement of the MCS in G2
    for k in range(c2):
        # get node from G2
        n2 = comp_mcs_G2[k]
        # copy alignment column
        sub_alignment_columns[current_alignment_nodes["cG2"][n2]] = deepcopy(columns_G2[n2])

    # end of function
    return(sub_alignment_columns)





# functions - ambiguous neighbors ##############################################





# function: --------------------------------------------------------------------
def progressive_alignment_get_ambiguous_neighbors(alignment_graph, alignment_columns, aligned_leaves):

    # local variables (cython)
    cdef int i = 0
    cdef int k = 0
    cdef int entry = 0

    # local variables (python)
    ambiguous = True
    cdef list nodes = []
    cdef tuple pair = tuple()
    cdef dict alignment_ambiguous_neighbors = dict()

    # get vertices in alignment
    nodes = deepcopy(list(alignment_graph.nodes(data = False)))

    # get ambiguous pairs
    for pair in combinations(nodes, r = 2):

        # continue only if not neighbors
        if((not alignment_graph.has_edge(pair[0], pair[1])) and (not alignment_graph.has_edge(pair[1], pair[0]))):

            # reinitialize variables
            ambiguous = True

            # compare column entries
            for i in range(len(aligned_leaves)):
                # get valid entry
                entry = aligned_leaves[i]
                # test disjointness in valid entry
                if((not alignment_columns[pair[0]][entry] == None) and (not alignment_columns[pair[1]][entry] == None)):
                    ambiguous = False
                    break

            # save if ambiguous neighbors
            if(ambiguous):

                # first node
                if(pair[0] in list(alignment_ambiguous_neighbors.keys())):
                    alignment_ambiguous_neighbors[pair[0]].append(deepcopy(pair[1]))
                else:
                    alignment_ambiguous_neighbors[pair[0]] = [deepcopy(pair[1])]

                # second node
                if(pair[1] in list(alignment_ambiguous_neighbors.keys())):
                    alignment_ambiguous_neighbors[pair[1]].append(deepcopy(pair[0]))
                else:
                    alignment_ambiguous_neighbors[pair[1]] = [deepcopy(pair[0])]

    # end of function
    return(alignment_ambiguous_neighbors)





# functions - consensus graphs #################################################





# function: --------------------------------------------------------------------
def progressive_alignment_consensus_graphs(alignment_graph, alignment_columns, N):

    # local variables (cython)
    cdef int p = 0
    cdef int k = 0
    cdef int order = 0
    cdef int gap_count = 0

    # local variables (python)
    cdef list bad_nodes = []
    cdef list alignment_nodes = []
    cdef dict not_gap_count = dict()
    cdef dict all_consensus_graphs = dict()
    consensus_graph = None

    # get vertices of alignment
    alignment_nodes = deepcopy(list(alignment_graph.nodes()))
    order = len(alignment_nodes)

    # get presence of vertices of alignment in input graphs
    for k in range(order):
        gap_count = alignment_columns[alignment_nodes[k]].count(None)
        not_gap_count[alignment_nodes[k]] = N - gap_count

    # get the consensus graph corresponding to each presence value
    for p in range(1, N+1):

        # reinitialize containers
        bad_nodes = []

        # get vertices failing to have not_gap_count at least p
        for k in range(N):
            if(not_gap_count[alignment_nodes[k]] < p):
                bad_nodes.append(alignment_nodes[k])

        # induce graph from vertices having not_gap_count >= p
        consensus_graph = deepcopy(alignment_graph)
        consensus_graph.remove_nodes_from(bad_nodes)
        all_consensus_graphs[str(p) + "/" + str(N)] = deepcopy(consensus_graph)

    # end of function
    return(all_consensus_graphs)





################################################################################
################################################################################
