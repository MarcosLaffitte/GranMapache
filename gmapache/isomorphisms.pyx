################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: isomorphisms                                                       #
#                                                                              #
# - Description: analysis of isomorphisms and automorphisms.                   #
#                                                                              #
# - NOTES:                                                                     #
#                                                                              #
#   * here we refer as "ring" to what in the literature of the VF2 is known    #
#     instead as "terminal" sets, for the simple reason that these sets look   #
#     like rings/doughnuts around the set of matched vertices when drawing     #
#     them on the whiteboard in a concentric layout.                           #
#                                                                              #
#   * the set of neighbors of a vertex that are not in the match and not in    #
#     the ring is called the set of "extern" neighbors.                        #
#                                                                              #
#   * the data structures holding the matched vertices, the rings and also     #
#     the unmatched vertices have both unordered and ordered representations   #
#     so as to make the most out of (average) constant look-ups with unordered #
#     structures, contiguous memory traversals with vectors, and as well as    #
#     for some linked-lists advantages.                                        #
#                                                                              #
#   * inside the intensive routines and in associated search-parameters, the   #
#     domain graph is called G, while the codomain graph is always called H.   #
#                                                                              #
################################################################################





# C++ string encoding ##########################################################
# cython: c_string_type=unicode, c_string_encoding=utf8





# dependencies #################################################################





# already in python ------------------------------------------------------------
from copy import deepcopy
from sys import getrecursionlimit, setrecursionlimit





# not in python ----------------------------------------------------------------
import networkx as nx





# cython specifics -------------------------------------------------------------
import cython
from libcpp cimport bool as cpp_bool
from libcpp.set cimport set as cpp_set
from libcpp.pair cimport pair as cpp_pair
from libcpp.list cimport list as cpp_list
from libcpp.stack cimport stack as cpp_stack
from libcpp.vector cimport vector as cpp_vector
from libcpp.string cimport string as cpp_string
from libcpp.unordered_set cimport unordered_set as cpp_unordered_set
from libcpp.unordered_map cimport unordered_map as cpp_unordered_map
cdef extern from "<algorithm>" namespace "std":
    # turn integer to string
    cpp_string to_string(int value)
    # get the ceiling of a number
    float  ceil(float num)
    # sort a vector
    void sort[Iter](Iter first, Iter last)
    # reverse vector
    void reverse[Iter](Iter first, Iter last)
    # find element in vector
    Iter find[Iter, Const](Iter first, Iter last, Const value)





# custom dependencies ----------------------------------------------------------
from .integerization import encode_graphs, decode_graphs, decode_match





# C/C++ structs ################################################################





# structs - graphs #############################################################





# struct: undirected graph -----------------------------------------------------
cdef struct isomorphisms_undirected_graph:
    # node data
    cpp_unordered_map[int, int] nodes
    # loop data
    cpp_unordered_set[int] loops
    # edge data
    cpp_unordered_map[cpp_string, int] edges
    # neighbors data (for constant look-ups)
    cpp_unordered_map[int, cpp_unordered_set[int]] neighbors
    # neighbors data complement (for constant look-ups)
    cpp_unordered_map[int, cpp_unordered_set[int]] neighbors_complement
    # ordered neighbors data (for insertion into ring, from original or complement)
    cpp_unordered_map[int, cpp_vector[int]] neighbors_ordered





# struct: directed graph -------------------------------------------------------
cdef struct isomorphisms_directed_graph:
    # node data
    cpp_unordered_map[int, int] nodes
    # loop data
    cpp_unordered_set[int] loops
    # edge data
    cpp_unordered_map[cpp_string, int] edges
    # neighbors data (for constant look-ups)
    cpp_unordered_map[int, cpp_unordered_set[int]] in_neighbors
    cpp_unordered_map[int, cpp_unordered_set[int]] out_neighbors
    # neighbors data complement (for constant look-ups)
    cpp_unordered_map[int, cpp_unordered_set[int]] in_neighbors_complement
    cpp_unordered_map[int, cpp_unordered_set[int]] out_neighbors_complement
    # ordered neighbors data (for insertion into ring, from original or complement)
    cpp_unordered_map[int, cpp_vector[int]] in_neighbors_ordered
    cpp_unordered_map[int, cpp_vector[int]] out_neighbors_ordered
    # neighbors in underlying undirected graph
    cpp_unordered_map[int, cpp_unordered_set[int]] connectivity_neighbors





# structs - states of search space  ############################################





# structure: state for the VF2-like isomorphism search - undirected ------------
cdef struct isomorphisms_state_undirected:
    # match information
    cpp_set[cpp_pair[int, int]] match
    cpp_unordered_set[int] match_G
    cpp_unordered_set[int] match_H
    cpp_unordered_map[int, int] forward_match
    cpp_unordered_map[int, int] inverse_match
    # ring information
    cpp_unordered_set[int] ring_G
    cpp_unordered_set[int] ring_H
    # extern information
    cpp_unordered_set[int] unmatched_G
    cpp_unordered_set[int] unmatched_H
    # ordered candidates
    cpp_list[int] ring_G_ordered
    cpp_list[int] ring_H_ordered
    cpp_list[int] unmatched_G_ordered
    cpp_list[int] unmatched_H_ordered





# structure: state for the VF2-like isomorphism search - directed --------------
cdef struct isomorphisms_state_directed:
    # match information
    cpp_set[cpp_pair[int, int]] match
    cpp_unordered_set[int] match_G
    cpp_unordered_set[int] match_H
    cpp_unordered_map[int, int] forward_match
    cpp_unordered_map[int, int] inverse_match
    # ring information
    cpp_unordered_set[int] in_ring_G
    cpp_unordered_set[int] in_ring_H
    cpp_unordered_set[int] out_ring_G
    cpp_unordered_set[int] out_ring_H
    # extern information
    cpp_unordered_set[int] unmatched_G
    cpp_unordered_set[int] unmatched_H
    # ordered candidates
    cpp_list[int] in_ring_G_ordered
    cpp_list[int] in_ring_H_ordered
    cpp_list[int] out_ring_G_ordered
    cpp_list[int] out_ring_H_ordered
    cpp_list[int] unmatched_G_ordered
    cpp_list[int] unmatched_H_ordered





# structs - parameters of search ###############################################





# struct: overall search parameters and constant information -------------------
cdef struct isomorphisms_search_params:
    # use node labels
    cpp_bool node_labels
    # use edge labels
    cpp_bool edge_labels
    # directed or undirected graphs
    cpp_bool directed_graphs
    # return all isomorphisms
    cpp_bool all_isomorphisms
    # analyzing complement
    cpp_bool complement
    # expected amount of matches (order of the graphs)
    size_t expected_order
    int expected_order_int
    # total order for VF2-like search
    cpp_unordered_map[int, int] total_order_G
    cpp_unordered_map[int, int] total_order_H
    # inverse total order for VF2-like search
    cpp_unordered_map[int, int] inverse_total_order_G
    cpp_unordered_map[int, int] inverse_total_order_H
    # information of anchor induced by degrees and / or labels
    cpp_unordered_set[int] induced_anchor_G
    cpp_unordered_set[int] induced_anchor_H
    cpp_vector[cpp_pair[int, int]] induced_encoded_anchor





# structs - auxiliary structs for search #######################################





# struct: nodes added or removed from rings during recursive search ------------
cdef struct isomorphisms_change_in_state_undirected:
    # node effectively removed from ring in G
    cpp_bool removed_node_ring_G
    # node effectively removed from ring in H
    cpp_bool removed_node_ring_H
    # neighbors effectively added to ring in G
    cpp_stack[int] added_neighbors_ring_G
    # neighbors effectively added to ring in H
    cpp_stack[int] added_neighbors_ring_H





# struct: nodes added or removed from rings during recursive search ------------
cdef struct isomorphisms_change_in_state_directed:
    # node effectively removed from ring in G
    cpp_bool removed_node_in_ring_G
    cpp_bool removed_node_out_ring_G
    # node effectively removed from ring in H
    cpp_bool removed_node_in_ring_H
    cpp_bool removed_node_out_ring_H
    # neighbors effectively added to ring in G
    cpp_stack[int] added_neighbors_in_ring_G
    cpp_stack[int] added_neighbors_out_ring_G
    # neighbors effectively added to ring in H
    cpp_stack[int] added_neighbors_in_ring_H
    cpp_stack[int] added_neighbors_out_ring_H





# algorithms ###################################################################





# functions - search isomorphisms - wrapper ####################################





# function: callable wrapper for searching for isomorphisms --------------------
def search_isomorphisms(nx_G = nx.Graph(),           # can also be a networkx DiGraph
                        nx_H = nx.Graph(),           # can also be a networkx DiGraph
                        node_labels = False,         # consider node labels when evaluating isomorphisms
                        edge_labels = False,         # consider edge labels when evaluating isomorphisms
                        all_isomorphisms = False):   # by default stops when finding one isomorphism (if any)

    # description
    """
    > description:
    receives two networkx (di-)graphs G and H both directed or both undirected,
    possibly but not necessarily with either node labels and/or edge labels, with
    the same number of nodes and edges, and a boolean variable indicating if the
    function should finish when finding only one isomorphism from G to H (if any),
    or if it should search for all possible isomorphisms from G to H and return
    them. In addition, the VF2-like search requires a total order for the nodes
    of one of the input graphs (traditionally the domain graph). In principle this
    can be an arbitrary total order, though it has been shown that certain orders
    may improve the search (see VF2++). Here we order the nodes in both graphs
    based on their (out) degree in descending order by default.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * node_labels - boolean indicating if all node labels should be considered for
    the search or if they should be ignored (default).
    * edge_labels - boolean indicating if all edge labels should be considered for
    the search or if they should be ignored (default).
    * all_isomorphisms - boolean variable indicating if the function should stop
    as soon as one isomorphism is found (if any) -default behavior- or if it
    should search for all possible isomorphisms between the graphs.

    > output:
    * isomorphisms - (possibly empty) list of isomorphisms, each as a list of
    2-tuples (x, y) of nodes x from G and y from H representing the one-to-one
    correspondences preserving adjacency and labels.
    * found_isomorphism - boolean value indicating if any isomorphism was found
    and equivalently if G and H are isomorphic (labeled) (di-)graphs.

    > calls:
    * gmapache.integerization.encode_graphs
    * gmapache.integerization.decode_graphs
    * gmapache.integerization.decode_match
    * gmapache.isomorphisms.search_isomorphisms_input_correctness
    * gmapache.isomorphisms.search_isomorphisms_order_nodes_by_degree
    * gmapache.isomorphisms.search_isomorphisms_label_consistency
    * gmapache.isomorphisms.search_isomorphisms_order_neighbors
    * gmapache.isomorphisms.search_isomorphisms_undirected
    * gmapache.isomorphisms.search_isomorphisms_directed
    """

    # output holders
    found_isomorphism = False
    cdef list isomorphisms = []

    # fixed threshold parameters
    cdef float limit_edges = 0.98
    cdef float bigger_graphs = 15
    cdef float scalation_value = 1.5

    # local variables (cython)
    cdef int i = 0
    cdef int deg = 0
    cdef int node = 0
    cdef int deg1 = 0
    cdef int deg2 = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int label = 0
    cdef int counter = 0
    cdef int current_limit = 0
    cdef int required_limit = 0
    cdef cpp_bool reachable = True
    cdef cpp_bool input_correctness = True
    cdef cpp_bool consistent_labels = True
    cdef cpp_bool feasible_induced_anchor = True
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[cpp_pair[int, int], cpp_pair[int, int]] each_info_tuple
    cdef cpp_vector[int] deg_G
    cdef cpp_vector[int] deg_H
    cdef cpp_vector[int] in_deg_G
    cdef cpp_vector[int] in_deg_H
    cdef cpp_vector[int] out_deg_G
    cdef cpp_vector[int] out_deg_H
    cdef cpp_vector[int] each_vector
    cdef cpp_vector[cpp_vector[int]] all_edges_G
    cdef cpp_vector[cpp_vector[int]] all_edges_H
    cdef cpp_vector[cpp_pair[int, int]] all_nodes_G
    cdef cpp_vector[cpp_pair[int, int]] all_nodes_H
    cdef cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] temp_nodes_G
    cdef cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] temp_nodes_H
    cdef cpp_vector[cpp_set[cpp_pair[int, int]]] encoded_isomorphisms
    cdef cpp_set[cpp_pair[int, int]] each_isomorphism
    cdef cpp_unordered_map[int, int] node_degrees
    cdef isomorphisms_search_params params
    cdef isomorphisms_directed_graph directed_G
    cdef isomorphisms_directed_graph directed_H
    cdef isomorphisms_undirected_graph undirected_G
    cdef isomorphisms_undirected_graph undirected_H
    cdef isomorphisms_state_directed initial_state_directed
    cdef isomorphisms_state_undirected initial_state_undirected

    # local variables (python)
    cdef list encoded_graphs = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_labels = dict()
    cdef dict encoded_edge_labels = dict()
    node_obj = None
    complement_G = None
    complement_H = None
    undirected_copy = None

    # test input correctness
    input_correctness = search_isomorphisms_input_correctness(nx_G, nx_H, node_labels, edge_labels, all_isomorphisms)
    if(not input_correctness):
        return([], False)

    # quick test by comparing the degree sequences of the input graphs
    params.directed_graphs = nx.is_directed(nx_G)
    if(params.directed_graphs):

        # get in-degrees
        in_deg_G = [deg for (node_obj, deg) in list(nx_G.in_degree())]
        in_deg_H = [deg for (node_obj, deg) in list(nx_H.in_degree())]

        # sort in-degrees
        sort(in_deg_G.begin(), in_deg_G.end())
        sort(in_deg_H.begin(), in_deg_H.end())

        # compare in-degree sequences
        if(not in_deg_G == in_deg_H):
            return([], False)

        # get out-degrees
        out_deg_G = [deg for (node_obj, deg) in list(nx_G.out_degree())]
        out_deg_H = [deg for (node_obj, deg) in list(nx_H.out_degree())]

        # sort out-degrees
        sort(out_deg_G.begin(), out_deg_G.end())
        sort(out_deg_H.begin(), out_deg_H.end())

        # compare out-degree sequences
        if(not out_deg_G == out_deg_H):
            return([], False)

    else:

        # get degrees
        deg_G = [deg for (node_obj, deg) in list(nx_G.degree())]
        deg_H = [deg for (node_obj, deg) in list(nx_H.degree())]

        # sort degrees
        sort(deg_G.begin(), deg_G.end())
        sort(deg_H.begin(), deg_H.end())

        # compare degree sequences
        if(not deg_G == deg_H):
            return([], False)

    # save input parameters
    params.node_labels = node_labels
    params.edge_labels = edge_labels
    params.all_isomorphisms = all_isomorphisms
    params.expected_order = nx_H.order()
    params.expected_order_int = nx_H.order()
    if(params.directed_graphs):
        params.complement = (nx_G.size() > (ceil((params.expected_order * (params.expected_order-1)) * limit_edges) + params.expected_order)) and (params.expected_order >= bigger_graphs)
    else:
        params.complement = (nx_G.size() > (ceil((params.expected_order * (params.expected_order-1)/2) * limit_edges) + params.expected_order)) and (params.expected_order >= bigger_graphs)

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_labels, encoded_edge_labels = encode_graphs([nx_G, nx_H])

    # get complement if necessary
    if(params.complement):
        complement_G = nx.complement(encoded_graphs[0])
        complement_H = nx.complement(encoded_graphs[1])

    # prototype list of nodes
    if(params.directed_graphs):

        # get node information for directed G
        if(params.complement):
            temp_nodes_G = [((node, info["GMNL"]), (complement_G.out_degree[node], complement_G.in_degree[node])) for (node, info) in encoded_graphs[0].nodes(data = True)]
        else:
            temp_nodes_G = [((node, info["GMNL"]), (encoded_graphs[0].out_degree[node], encoded_graphs[0].in_degree[node])) for (node, info) in encoded_graphs[0].nodes(data = True)]

        # get node information for directed H
        if(params.complement):
            temp_nodes_H = [((node, info["GMNL"]), (complement_H.out_degree[node], complement_H.in_degree[node])) for (node, info) in encoded_graphs[1].nodes(data = True)]
        else:
            temp_nodes_H = [((node, info["GMNL"]), (encoded_graphs[1].out_degree[node], encoded_graphs[1].in_degree[node])) for (node, info) in encoded_graphs[1].nodes(data = True)]

    else:

        # get node information for undirected G
        if(params.complement):
            temp_nodes_G = [((node, info["GMNL"]), (complement_G.degree[node], 0)) for (node, info) in encoded_graphs[0].nodes(data = True)]
        else:
            temp_nodes_G = [((node, info["GMNL"]), (encoded_graphs[0].degree[node], 0)) for (node, info) in encoded_graphs[0].nodes(data = True)]

        # get node information for undirected H
        if(params.complement):
            temp_nodes_H = [((node, info["GMNL"]), (complement_H.degree[node], 0)) for (node, info) in encoded_graphs[1].nodes(data = True)]
        else:
            temp_nodes_H = [((node, info["GMNL"]), (encoded_graphs[1].degree[node], 0)) for (node, info) in encoded_graphs[1].nodes(data = True)]

    # prepare nodes and neighbors
    if(params.directed_graphs):

        # nodes of directed G
        directed_G.loops = list(nx.nodes_with_selfloops(encoded_graphs[0]))
        for each_info_tuple in temp_nodes_G:
            # save node
            directed_G.nodes.insert(each_info_tuple.first)
            # neighbors for G
            if((not params.complement) or params.node_labels or params.edge_labels):
                directed_G.in_neighbors[each_info_tuple.first.first] = set(encoded_graphs[0].predecessors(each_info_tuple.first.first))
                directed_G.out_neighbors[each_info_tuple.first.first] = set(encoded_graphs[0].neighbors(each_info_tuple.first.first))
            # neighbors for complement of G
            if(params.complement):
                directed_G.in_neighbors_complement[each_info_tuple.first.first] = set(complement_G.predecessors(each_info_tuple.first.first))
                directed_G.out_neighbors_complement[each_info_tuple.first.first] = set(complement_G.neighbors(each_info_tuple.first.first))

        # nodes of directed H
        directed_H.loops = list(nx.nodes_with_selfloops(encoded_graphs[1]))
        for each_info_tuple in temp_nodes_H:
            # save node
            directed_H.nodes.insert(each_info_tuple.first)
            # neighbors for H
            if((not params.complement) or params.node_labels or params.edge_labels):
                directed_H.in_neighbors[each_info_tuple.first.first] = set(encoded_graphs[1].predecessors(each_info_tuple.first.first))
                directed_H.out_neighbors[each_info_tuple.first.first] = set(encoded_graphs[1].neighbors(each_info_tuple.first.first))
            # neighbors for complement of H
            if(params.complement):
                directed_H.in_neighbors_complement[each_info_tuple.first.first] = set(complement_H.predecessors(each_info_tuple.first.first))
                directed_H.out_neighbors_complement[each_info_tuple.first.first] = set(complement_H.neighbors(each_info_tuple.first.first))

    else:

        # nodes of undirected G
        undirected_G.loops = list(nx.nodes_with_selfloops(encoded_graphs[0]))
        for each_info_tuple in temp_nodes_G:
            # save node
            undirected_G.nodes.insert(each_info_tuple.first)
            # neighbors for G
            if((not params.complement) or params.node_labels or params.edge_labels):
                undirected_G.neighbors[each_info_tuple.first.first] = set(encoded_graphs[0].neighbors(each_info_tuple.first.first))
            # neighbors for complement of G
            if(params.complement):
                undirected_G.neighbors_complement[each_info_tuple.first.first] = set(complement_G.neighbors(each_info_tuple.first.first))

        # nodes of undirected H
        undirected_H.loops = list(nx.nodes_with_selfloops(encoded_graphs[1]))
        for each_info_tuple in temp_nodes_H:
            # save node
            undirected_H.nodes.insert(each_info_tuple.first)
            # neighbors for H
            if((not params.complement) or params.node_labels or params.edge_labels):
                undirected_H.neighbors[each_info_tuple.first.first] = set(encoded_graphs[1].neighbors(each_info_tuple.first.first))
            # neighbors for complement of H
            if(params.complement):
                undirected_H.neighbors_complement[each_info_tuple.first.first] = set(complement_H.neighbors(each_info_tuple.first.first))

    # prepare edges
    if(params.edge_labels):
        if(params.directed_graphs):

            # edges for G
            all_edges_G = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[0].edges(data = True)]
            for each_vector in all_edges_G:
                temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                directed_G.edges[temp_str] = each_vector[2]

            # edges for H
            all_edges_H = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[1].edges(data = True)]
            for each_vector in all_edges_H:
                temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                directed_H.edges[temp_str] = each_vector[2]

        else:

            # edges for G
            all_edges_G = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[0].edges(data = True)]
            for each_vector in all_edges_G:
                if(each_vector[0] == each_vector[1]):
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[0])
                    undirected_G.edges[temp_str] = each_vector[2]
                else:
                    # save the two label edges to simplify future access
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                    undirected_G.edges[temp_str] = each_vector[2]
                    temp_str = to_string(each_vector[1]) + comma + to_string(each_vector[0])
                    undirected_G.edges[temp_str] = each_vector[2]

            # edges for H
            all_edges_H = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[1].edges(data = True)]
            for each_vector in all_edges_H:
                if(each_vector[0] == each_vector[1]):
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[0])
                    undirected_H.edges[temp_str] = each_vector[2]
                else:
                    # save the two label edges to simplify future access
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                    undirected_H.edges[temp_str] = each_vector[2]
                    temp_str = to_string(each_vector[1]) + comma + to_string(each_vector[0])
                    undirected_H.edges[temp_str] = each_vector[2]

    # test consistency of labels if required and determine label-induced anchor (if any)
    if(params.node_labels or params.edge_labels):

        # test consistency of node and/or edge labels
        if(params.directed_graphs):
            consistent_labels = search_isomorphisms_label_consistency(params.node_labels, params.edge_labels,
                                                                      directed_G.nodes, directed_H.nodes,
                                                                      directed_G.edges, directed_H.edges,
                                                                      params)
        else:
            consistent_labels = search_isomorphisms_label_consistency(params.node_labels, params.edge_labels,
                                                                      undirected_G.nodes, undirected_H.nodes,
                                                                      undirected_G.edges, undirected_H.edges,
                                                                      params)

        # if labels are not consistent then the graphs cannot be isomorphic
        if(not consistent_labels):
            return([], False)

    # determine degree-induced anchor (if any and only if no label-check is required)
    if((not params.node_labels) and (params.induced_encoded_anchor.empty())):
        feasible_induced_anchor = search_isomorphisms_test_degree_induced_anchor(temp_nodes_G, temp_nodes_H, params)
        # for directed graphs, if degrees induce an inconsistent anchor then the graphs cannot be isomorphic
        if(not feasible_induced_anchor):
            return([], False)

    # determine and create total orders
    if(not params.induced_encoded_anchor.empty()):
        if(params.directed_graphs):

            # get concentric order for directed G
            if(params.complement):
                undirected_copy = deepcopy(complement_G)
            else:
                undirected_copy = deepcopy(encoded_graphs[0])
            undirected_copy = undirected_copy.to_undirected()
            directed_G.connectivity_neighbors = {node:set(undirected_copy.neighbors(node)) for ((node, label), (deg1, deg2)) in temp_nodes_G}
            node_degrees = {node:deg1 for ((node, label), (deg1, deg2)) in temp_nodes_G}
            reachable = search_isomorphisms_order_nodes_by_concentric_reachability(params.induced_anchor_G, temp_nodes_G, directed_G.connectivity_neighbors, node_degrees, all_nodes_G, directed_G.nodes)

            # get concentric order for directed H
            if(params.complement):
                undirected_copy = deepcopy(complement_H)
            else:
                undirected_copy = deepcopy(encoded_graphs[1])
            undirected_copy = undirected_copy.to_undirected()
            directed_H.connectivity_neighbors = {node:set(undirected_copy.neighbors(node)) for ((node, label), (deg1, deg2)) in temp_nodes_H}
            node_degrees = {node:deg1 for ((node, label), (deg1, deg2)) in temp_nodes_H}
            reachable = search_isomorphisms_order_nodes_by_concentric_reachability(params.induced_anchor_H, temp_nodes_H, directed_H.connectivity_neighbors, node_degrees, all_nodes_H, directed_H.nodes)

        else:

            # get concentric order for undirected G
            node_degrees = {node:deg1 for ((node, label), (deg1, deg2)) in temp_nodes_G}
            if(params.complement):
                reachable = search_isomorphisms_order_nodes_by_concentric_reachability(params.induced_anchor_G, temp_nodes_G, undirected_G.neighbors_complement, node_degrees, all_nodes_G, undirected_G.nodes)
            else:
                reachable = search_isomorphisms_order_nodes_by_concentric_reachability(params.induced_anchor_G, temp_nodes_G, undirected_G.neighbors, node_degrees, all_nodes_G, undirected_G.nodes)

            # get concentric order for undirected H
            node_degrees = {node:deg1 for ((node, label), (deg1, deg2)) in temp_nodes_H}
            if(params.complement):
                reachable = search_isomorphisms_order_nodes_by_concentric_reachability(params.induced_anchor_H, temp_nodes_H, undirected_H.neighbors_complement, node_degrees, all_nodes_H, undirected_H.nodes)
            else:
                reachable = search_isomorphisms_order_nodes_by_concentric_reachability(params.induced_anchor_H, temp_nodes_H, undirected_H.neighbors, node_degrees, all_nodes_H, undirected_H.nodes)

    else:

        # get degree-based order for directed or undirected G
        search_isomorphisms_order_nodes_by_degree(temp_nodes_G, all_nodes_G)
        # get degree-based order for directed or undirected H
        search_isomorphisms_order_nodes_by_degree(temp_nodes_H, all_nodes_H)

    # store total orders
    if(params.directed_graphs):

        # total order of directed G
        counter = 0
        for each_pair in all_nodes_G:
            counter = counter + 1
            params.total_order_G[each_pair.first] = counter
            params.inverse_total_order_G[counter] = each_pair.first
            initial_state_directed.unmatched_G.insert(each_pair.first)
            initial_state_directed.unmatched_G_ordered.push_back(each_pair.first)

        # total order of directed H
        counter = 0
        for each_pair in all_nodes_H:
            counter = counter + 1
            params.total_order_H[each_pair.first] = counter
            params.inverse_total_order_H[counter] = each_pair.first
            initial_state_directed.unmatched_H.insert(each_pair.first)
            initial_state_directed.unmatched_H_ordered.push_back(each_pair.first)

    else:

        # total order of undirected G
        counter = 0
        for each_pair in all_nodes_G:
            counter = counter + 1
            params.total_order_G[each_pair.first] = counter
            params.inverse_total_order_G[counter] = each_pair.first
            initial_state_undirected.unmatched_G.insert(each_pair.first)
            initial_state_undirected.unmatched_G_ordered.push_back(each_pair.first)

        # total order of undirected H
        counter = 0
        for each_pair in all_nodes_H:
            counter = counter + 1
            params.total_order_H[each_pair.first] = counter
            params.inverse_total_order_H[counter] = each_pair.first
            initial_state_undirected.unmatched_H.insert(each_pair.first)
            initial_state_undirected.unmatched_H_ordered.push_back(each_pair.first)

    # prepare ordered neighbors
    if(params.complement):

        if(params.directed_graphs):

            # order nodes of G
            search_isomorphisms_order_neighbors(all_nodes_G, directed_G.in_neighbors_complement, directed_G.in_neighbors_ordered,
                                                params.total_order_G, params.inverse_total_order_G)
            search_isomorphisms_order_neighbors(all_nodes_G, directed_G.out_neighbors_complement, directed_G.out_neighbors_ordered,
                                                params.total_order_G, params.inverse_total_order_G)

            # order nodes of H
            search_isomorphisms_order_neighbors(all_nodes_H, directed_H.in_neighbors_complement, directed_H.in_neighbors_ordered,
                                                params.total_order_H, params.inverse_total_order_H)
            search_isomorphisms_order_neighbors(all_nodes_H, directed_H.out_neighbors_complement, directed_H.out_neighbors_ordered,
                                                params.total_order_H, params.inverse_total_order_H)

        else:

            # order nodes of G
            search_isomorphisms_order_neighbors(all_nodes_G, undirected_G.neighbors_complement, undirected_G.neighbors_ordered,
                                                params.total_order_G, params.inverse_total_order_G)

            # order nodes of H
            search_isomorphisms_order_neighbors(all_nodes_H, undirected_H.neighbors_complement, undirected_H.neighbors_ordered,
                                                params.total_order_H, params.inverse_total_order_H)

    else:

        if(params.directed_graphs):

            # order nodes of G
            search_isomorphisms_order_neighbors(all_nodes_G, directed_G.in_neighbors, directed_G.in_neighbors_ordered,
                                                params.total_order_G, params.inverse_total_order_G)
            search_isomorphisms_order_neighbors(all_nodes_G, directed_G.out_neighbors, directed_G.out_neighbors_ordered,
                                                params.total_order_G, params.inverse_total_order_G)

            # order nodes of H
            search_isomorphisms_order_neighbors(all_nodes_H, directed_H.in_neighbors, directed_H.in_neighbors_ordered,
                                                params.total_order_H, params.inverse_total_order_H)
            search_isomorphisms_order_neighbors(all_nodes_H, directed_H.out_neighbors, directed_H.out_neighbors_ordered,
                                                params.total_order_H, params.inverse_total_order_H)

        else:

            # order nodes of G
            search_isomorphisms_order_neighbors(all_nodes_G, undirected_G.neighbors, undirected_G.neighbors_ordered,
                                                params.total_order_G, params.inverse_total_order_G)

            # order nodes of H
            search_isomorphisms_order_neighbors(all_nodes_H, undirected_H.neighbors, undirected_H.neighbors_ordered,
                                                params.total_order_H, params.inverse_total_order_H)

    # set recursion limit
    required_limit = nx_G.order()
    current_limit = getrecursionlimit()
    if(current_limit < (scalation_value * required_limit)):
        setrecursionlimit(int(scalation_value * required_limit))

    # prepare initial state with induced anchor (if any)
    if(not params.induced_encoded_anchor.empty()):
        if(params.directed_graphs):
            feasible_induced_anchor = search_isomorphisms_prepare_feasible_initial_state_directed(params, directed_G, directed_H, initial_state_directed)
        else:
            feasible_induced_anchor = search_isomorphisms_prepare_feasible_initial_state_undirected(params, undirected_G, undirected_H, initial_state_undirected)

    # evaluate isomorphisms
    if(feasible_induced_anchor):
        if(params.directed_graphs):
            search_isomorphisms_directed(params, initial_state_directed, directed_G, directed_H, encoded_isomorphisms)
        else:
            search_isomorphisms_undirected(params, initial_state_undirected, undirected_G, undirected_H, encoded_isomorphisms)

    # decode isomorphisms
    for each_isomorphism in encoded_isomorphisms:
        isomorphisms.append(decode_match(list(each_isomorphism), encoded_node_names))

    # check if the graphs were isomorphic
    if(len(isomorphisms) > 0):
        found_isomorphism = True

    # end of function
    return(isomorphisms, found_isomorphism)





# functions - input consistency ################################################





# function: test input correctness ---------------------------------------------
cdef cpp_bool search_isomorphisms_input_correctness(nx_G, nx_H, node_labels, edge_labels, all_isomorphisms):

    # local variables (cython)
    cdef int minimum_order = 0
    cdef int maximum_order = 0
    cdef cpp_bool inserted = False
    cdef cpp_unordered_set[int] all_orders

    # local variables (python)
    test_int = 0
    test_bool = False
    cdef list nodes_G = []
    cdef list nodes_H = []
    cdef list all_keys = []
    cdef dict test_dict = dict()
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    each_node = None

    # check that first argument is networkx graph or digraph
    if(type(nx_G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))

    # check that second argument is networkx graph or digraph
    if(type(nx_H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))

    # check that third argument is boolean
    if(type(node_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument node_labels must be a boolean variable."))

    # check that fourth argument is boolean
    if(type(edge_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument edge_labels must be a boolean variable."))

    # check that fifth argument is boolean
    if(type(all_isomorphisms) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument all_isomorphisms must be a boolean variable."))

    # check that input graphs are of the same type
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))

    # check that input graphs are not null graphs
    if((nx_G.order() == 0) or (nx_H.order() == 0)):
        raise(ValueError("gmapache: input graphs must have at least one node each."))

    # check that input graphs have the same number of nodes
    if(not nx_G.order() == nx_H.order()):
        return(False)

    # check that input graphs have the same number of edges
    if(not nx_G.size() == nx_H.size()):
        return(False)

    # end of function
    return(True)





# function: consistency of node or edge labels and label-induced anchor --------
cdef cpp_bool search_isomorphisms_label_consistency(cpp_bool & node_labels,
                                                    cpp_bool & edge_labels,
                                                    cpp_unordered_map[int, int] & nodes_G,
                                                    cpp_unordered_map[int, int] & nodes_H,
                                                    cpp_unordered_map[cpp_string, int] & edges_G,
                                                    cpp_unordered_map[cpp_string, int] & edges_H,
                                                    isomorphisms_search_params & params) noexcept:

    # local variables
    cdef int label = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[int, int] node_and_label
    cdef cpp_pair[int, int] label_and_count
    cdef cpp_pair[cpp_string, int] edge_and_label
    cdef cpp_vector[int] empty_vector
    cdef cpp_unordered_map[int, int] count_labels_G
    cdef cpp_unordered_map[int, int] count_labels_H
    cdef cpp_unordered_map[int, cpp_vector[int]] labels_to_nodes_G
    cdef cpp_unordered_map[int, cpp_vector[int]] labels_to_nodes_H

    # consistency of node labels and label-induced anchor
    if(node_labels):

        # get node label count on G
        for node_and_label in nodes_G:
            # get node label
            label = node_and_label.second
            # count node label
            if(count_labels_G.find(label) != count_labels_G.end()):
                count_labels_G[label] = count_labels_G[label] + 1
                labels_to_nodes_G[label].push_back(node_and_label.first)
            else:
                count_labels_G[label] = 1
                labels_to_nodes_G[label] = empty_vector
                labels_to_nodes_G[label].push_back(node_and_label.first)

        # get node label count on H
        for node_and_label in nodes_H:
            # get node label
            label = node_and_label.second
            # count node label
            if(count_labels_H.find(label) != count_labels_H.end()):
                count_labels_H[label] = count_labels_H[label] + 1
                labels_to_nodes_H[label].push_back(node_and_label.first)
            else:
                count_labels_H[label] = 1
                labels_to_nodes_H[label] = empty_vector
                labels_to_nodes_H[label].push_back(node_and_label.first)

        # compare node counts in both graphs per label
        if(count_labels_G.size() != count_labels_H.size()):
            return(False)
        for label_and_count in count_labels_G:
            if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                return(False)
            else:
                if(label_and_count.second != count_labels_H[label_and_count.first]):
                    return(False)
                else:
                    if(label_and_count.second == 1):
                        temp_pair.first = labels_to_nodes_G[label_and_count.first][0]
                        temp_pair.second = labels_to_nodes_H[label_and_count.first][0]
                        params.induced_anchor_G.insert(temp_pair.first)
                        params.induced_anchor_H.insert(temp_pair.second)
                        params.induced_encoded_anchor.push_back(temp_pair)

        # if there is only one node label then turn off node-label checking, since it is the same label
        if(count_labels_G.size() == 1):
            node_labels = False

    # consistency of edge labels
    if(edge_labels):

        # clear counts
        count_labels_G.clear()
        count_labels_H.clear()

        # get edge label count on G
        for edge_and_label in edges_G:
            # get edge label
            label = edge_and_label.second
            # count edge label
            if(count_labels_G.find(label) != count_labels_G.end()):
                count_labels_G[label] = count_labels_G[label] + 1
            else:
                count_labels_G[label] = 1

        # get edge label count on H
        for edge_and_label in edges_H:
            # get edge label
            label = edge_and_label.second
            # count edge label
            if(count_labels_H.find(label) != count_labels_H.end()):
                count_labels_H[label] = count_labels_H[label] + 1
            else:
                count_labels_H[label] = 1

        # compare edge counts in both graphs per label
        if(count_labels_G.size() != count_labels_H.size()):
            return(False)
        for label_and_count in count_labels_G:
            if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                return(False)
            else:
                if(label_and_count.second != count_labels_H[label_and_count.first]):
                    return(False)

        # if there is only one edge label then turn off edge-label checking, since it is the same label
        if(count_labels_G.size() == 1):
            edge_labels = False

    # end of function
    return(True)





# functions - input preparation ################################################





# function: test existence of degree-induced anchor ----------------------------
cdef cpp_bool search_isomorphisms_test_degree_induced_anchor(cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] & nodes_info_G,
                                                             cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] & nodes_info_H,
                                                             isomorphisms_search_params & params) noexcept:

    # local variables
    cdef int deg = 0
    cdef int in_deg = 0
    cdef int out_deg = 0
    cdef int node = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_pair[cpp_pair[int, int], cpp_pair[int, int]] each_info_tuple
    cdef cpp_vector[int] all_degrees
    cdef cpp_vector[int] all_in_degrees
    cdef cpp_vector[int] all_out_degrees
    cdef cpp_vector[int] temp_vector
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_degree_G
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_in_degree_G
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_out_degree_G
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_degree_H
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_in_degree_H
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_out_degree_H

    # check degree for undirected graphs
    if(not params.directed_graphs):

        # save nodes in G by degree
        for each_info_tuple in nodes_info_G:

            # get degree
            node = each_info_tuple.first.first
            deg = each_info_tuple.second.first

            # determine if degree was already encountered
            if(nodes_by_degree_G.find(deg) != nodes_by_degree_G.end()):
                # save new pair associated to corresponding degree
                nodes_by_degree_G[deg].push_back(node)
            else:
                # save new degree and create vector of pairs with nodes
                all_degrees.push_back(deg)
                temp_vector.clear()
                temp_vector.push_back(node)
                nodes_by_degree_G[deg] = temp_vector

        # save nodes in H by degree
        for each_info_tuple in nodes_info_H:

            # get degree
            node = each_info_tuple.first.first
            deg = each_info_tuple.second.first

            # determine if degree was already encountered
            if(nodes_by_degree_H.find(deg) != nodes_by_degree_H.end()):
                # save new pair associated to corresponding degree
                nodes_by_degree_H[deg].push_back(node)
            else:
                # save new degree and create vector of pairs with nodes
                temp_vector.clear()
                temp_vector.push_back(node)
                nodes_by_degree_H[deg] = temp_vector

        # evaluate degree-induced anchor
        # NOTE: at this point the graphs have the same degree sequence
        for deg in all_degrees:
            if(nodes_by_degree_G[deg].size() == 1):
                temp_pair.first = nodes_by_degree_G[deg][0]
                temp_pair.second = nodes_by_degree_H[deg][0]
                params.induced_anchor_G.insert(temp_pair.first)
                params.induced_anchor_H.insert(temp_pair.second)
                params.induced_encoded_anchor.push_back(temp_pair)

    # check in and out degrees for directed graphs
    if(params.directed_graphs):

        # save nodes in G by in and out-degree
        for each_info_tuple in nodes_info_G:

            # get degrees
            node = each_info_tuple.first.first
            out_deg = each_info_tuple.second.first
            in_deg = each_info_tuple.second.second

            # determine if in-degree was already encountered
            if(nodes_by_in_degree_G.find(in_deg) != nodes_by_in_degree_G.end()):
                # save new pair associated to corresponding degree
                nodes_by_in_degree_G[in_deg].push_back(node)
            else:
                # save new degree and create vector of pairs with nodes
                temp_vector.clear()
                temp_vector.push_back(node)
                nodes_by_in_degree_G[in_deg] = temp_vector

            # determine if out-degree was already encountered
            if(nodes_by_out_degree_G.find(out_deg) != nodes_by_out_degree_G.end()):
                # save new pair associated to corresponding degree
                nodes_by_out_degree_G[out_deg].push_back(node)
            else:
                # save new degree and create vector of pairs with nodes
                temp_vector.clear()
                temp_vector.push_back(node)
                nodes_by_out_degree_G[out_deg] = temp_vector

        # save nodes in H by in and out-degree
        for each_info_tuple in nodes_info_H:

            # get degrees
            node = each_info_tuple.first.first
            out_deg = each_info_tuple.second.first
            in_deg = each_info_tuple.second.second

            # determine if in-degree was already encountered
            if(nodes_by_in_degree_H.find(in_deg) != nodes_by_in_degree_H.end()):
                # save new pair associated to corresponding degree
                nodes_by_in_degree_H[in_deg].push_back(node)
            else:
                # save new degree and create vector of pairs with nodes
                temp_vector.clear()
                temp_vector.push_back(node)
                nodes_by_in_degree_H[in_deg] = temp_vector

            # determine if out-degree was already encountered
            if(nodes_by_out_degree_H.find(out_deg) != nodes_by_out_degree_H.end()):
                # save new pair associated to corresponding degree
                nodes_by_out_degree_H[out_deg].push_back(node)
            else:
                # save new degree and create vector of pairs with nodes
                temp_vector.clear()
                temp_vector.push_back(node)
                nodes_by_out_degree_H[out_deg] = temp_vector

        # traverse nodes again testing uniqueness of degrees
        # NOTE: at this point the graphs have the same in-degree and out-degree sequence
        for each_info_tuple in nodes_info_G:

            # get degrees
            node = each_info_tuple.first.first
            out_deg = each_info_tuple.second.first
            in_deg = each_info_tuple.second.second

            # test uniqueness of degrees and match vertices if required
            # NOTE: we map via both type of degrees for syntactic feasibility
            if((nodes_by_in_degree_G[in_deg].size() == 1) and (nodes_by_out_degree_G[out_deg].size() == 1)):
                if(nodes_by_in_degree_H[in_deg][0] != nodes_by_out_degree_H[out_deg][0]):
                    return(False)
                else:
                    temp_pair.first = nodes_by_in_degree_G[in_deg][0]
                    temp_pair.second = nodes_by_in_degree_H[in_deg][0]
                    params.induced_anchor_G.insert(temp_pair.first)
                    params.induced_anchor_H.insert(temp_pair.second)
                    params.induced_encoded_anchor.push_back(temp_pair)

    # end of function
    return(True)





# function: order nodes by decreasing (out) degree -----------------------------
cdef void search_isomorphisms_order_nodes_by_degree(cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] & nodes_info,
                                                    cpp_vector[cpp_pair[int, int]] & all_nodes) noexcept:

    # local variables
    cdef int deg = 0
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[cpp_pair[int, int], cpp_pair[int, int]] each_info_tuple
    cdef cpp_vector[int] all_degrees
    cdef cpp_vector[cpp_pair[int, int]] temp_vector
    cdef cpp_unordered_map[int, cpp_vector[cpp_pair[int, int]]] nodes_by_degree

    # partition node information by node degree
    for each_info_tuple in nodes_info:
        # determine if degree was already encountered
        if(nodes_by_degree.find(each_info_tuple.second.first) != nodes_by_degree.end()):
            # save new pair associated to corresponding degree
            nodes_by_degree[each_info_tuple.second.first].push_back(each_info_tuple.first)
        else:
            # save new degree and create vector of pairs with node information
            all_degrees.push_back(each_info_tuple.second.first)
            temp_vector.clear()
            temp_vector.push_back(each_info_tuple.first)
            nodes_by_degree[each_info_tuple.second.first] = temp_vector

    # sort encountered degrees in decreasing order
    sort(all_degrees.begin(), all_degrees.end())
    reverse(all_degrees.begin(), all_degrees.end())

    # traverse nodes one to save them by decreasing degree
    for deg in all_degrees:
        for each_pair in nodes_by_degree[deg]:
            all_nodes.push_back(each_pair)

    # end of function





# function: order nodes in concentric sets around the anchor -------------------
cdef cpp_bool search_isomorphisms_order_nodes_by_concentric_reachability(cpp_unordered_set[int] & initial_level,
                                                                         cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] & temp_nodes,
                                                                         cpp_unordered_map[int, cpp_unordered_set[int]] & connectivity_neighbors,
                                                                         cpp_unordered_map[int, int] & node_degrees,
                                                                         cpp_vector[cpp_pair[int, int]] & all_nodes,
                                                                         cpp_unordered_map[int, int] node_labels) noexcept:

    # local variables
    cdef cpp_bool reachable = True
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[cpp_pair[int, int], cpp_pair[int, int]] each_info_tuple
    cdef cpp_pair[cpp_pair[int, int], cpp_pair[int, int]] temp_info_tuple
    cdef cpp_vector[cpp_pair[int, int]] next_level
    cdef cpp_vector[cpp_pair[int, int]] current_level
    cdef cpp_vector[cpp_pair[int, int]] unmatched_level
    cdef cpp_vector[cpp_pair[cpp_pair[int, int], cpp_pair[int, int]]] pre_level
    cdef cpp_unordered_map[int, cpp_bool] visited

    # prepare visited-nodes array
    for each_info_tuple in temp_nodes:
        node = each_info_tuple.first.first
        visited[node] = False

    # prepare nodes for current level
    for node in initial_level:
        temp_info_tuple.first.first = node
        temp_info_tuple.first.second = node_labels[node]
        temp_info_tuple.second.first = node_degrees[node]
        temp_info_tuple.second.second = 0
        pre_level.push_back(temp_info_tuple)
        visited[node] = True

    # order and save nodes in current level
    search_isomorphisms_order_nodes_by_degree(pre_level, current_level)
    for each_pair in current_level:
        all_nodes.push_back(each_pair)

    # multi-source DFS
    while(not current_level.empty()):
        # reinitialize next_level
        pre_level.clear()
        next_level.clear()
        # iterate getting immediate neighbors not already ordered
        for each_pair in current_level:
            node1 = each_pair.first
            for node2 in connectivity_neighbors[node1]:
                # only assign oreder of not visited yet
                if(not visited[node2]):
                    # prepare nodes for next level
                    temp_info_tuple.first.first = node2
                    temp_info_tuple.first.second = node_labels[node2]
                    temp_info_tuple.second.first = node_degrees[node2]
                    temp_info_tuple.second.second = 0
                    pre_level.push_back(temp_info_tuple)
                    visited[node2] = True
        # level management
        search_isomorphisms_order_nodes_by_degree(pre_level, next_level)
        for each_pair in next_level:
            all_nodes.push_back(each_pair)
        # update nodes to be ordered
        current_level.clear()
        current_level = next_level

    # test for unreachable nodes
    pre_level.clear()
    unmatched_level.clear()
    for each_info_tuple in temp_nodes:
        node = each_info_tuple.first.first
        if(not visited[node]):
            pre_level.push_back(each_info_tuple)
            visited[node] = True

    # enumerate unreachable nodes
    if(not pre_level.empty()):
        # set reachability to false
        reachable = False
        # order and save unmatched vertices
        search_isomorphisms_order_nodes_by_degree(pre_level, unmatched_level)
        for each_pair in unmatched_level:
            all_nodes.push_back(each_pair)

    # end of function
    return(reachable)





# function: order neighbors according to the total order -----------------------
cdef void search_isomorphisms_order_neighbors(cpp_vector[cpp_pair[int, int]] & all_nodes,
                                              cpp_unordered_map[int, cpp_unordered_set[int]] & all_neighbors,
                                              cpp_unordered_map[int, cpp_vector[int]] & all_neighbors_ordered,
                                              cpp_unordered_map[int, int] & total_order,
                                              cpp_unordered_map[int, int] & inverse_total_order) noexcept:

    # local variables
    cdef int each_order = 0
    cdef int each_neighbor = 0
    cdef cpp_vector[int] temp_vector
    cdef cpp_vector[int] empty_vector
    cdef cpp_pair[int, int] each_pair

    # order neighbors of each node
    for each_pair in all_nodes:
        # get total orders
        temp_vector.clear()
        for each_neighbor in all_neighbors[each_pair.first]:
            temp_vector.push_back(total_order[each_neighbor])
        # sort total orders increasingly
        sort(temp_vector.begin(), temp_vector.end())
        # get ordered neighbors
        all_neighbors_ordered[each_pair.first] = empty_vector
        for each_order in temp_vector:
            all_neighbors_ordered[each_pair.first].push_back(inverse_total_order[each_order])

    # end of function





# function: prepare initial state with induced anchor - undirected -------------
cdef cpp_bool search_isomorphisms_prepare_feasible_initial_state_undirected(isomorphisms_search_params & params,
                                                                            isomorphisms_undirected_graph & G,
                                                                            isomorphisms_undirected_graph & H,
                                                                            isomorphisms_state_undirected & current_state) noexcept:

    # local variables
    cdef cpp_bool inserted = False
    cdef cpp_bool feasible_induced_anchor = True
    cdef cpp_bool semantic_feasibility = True
    cdef cpp_bool syntactic_feasibility = True
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int taken_from = 0
    cdef cpp_pair[int, int] each_pair
    cdef isomorphisms_change_in_state_undirected change_in_state

    # extend initial search state with each pair in the anchor
    for each_pair in params.induced_encoded_anchor:

        # reinitialize feasibility
        semantic_feasibility = True
        syntactic_feasibility = True

        # split nodes
        node1 = each_pair.first
        node2 = each_pair.second

        # both vertices should be in the ring or outside of it
        if(current_state.ring_G.find(node1) != current_state.ring_G.end()):
            taken_from = 0
            if(current_state.ring_H.find(node2) == current_state.ring_H.end()):
                feasible_induced_anchor = False
                break
        else:
            taken_from = 1
            if(current_state.ring_H.find(node2) != current_state.ring_H.end()):
                feasible_induced_anchor = False
                break

        # evaluate syntactic feasibility (possibly over complement graph)
        if(params.complement):
            syntactic_feasibility = syntactic_feasibility_undirected(node1,
                                                                     node2,
                                                                     current_state.ring_G,
                                                                     current_state.ring_H,
                                                                     current_state.match_G,
                                                                     current_state.match_H,
                                                                     G.loops,
                                                                     H.loops,
                                                                     current_state.forward_match,
                                                                     current_state.inverse_match,
                                                                     G.neighbors_complement,
                                                                     H.neighbors_complement,
                                                                     G.neighbors_ordered,
                                                                     H.neighbors_ordered)
        else:
            syntactic_feasibility = syntactic_feasibility_undirected(node1,
                                                                     node2,
                                                                     current_state.ring_G,
                                                                     current_state.ring_H,
                                                                     current_state.match_G,
                                                                     current_state.match_H,
                                                                     G.loops,
                                                                     H.loops,
                                                                     current_state.forward_match,
                                                                     current_state.inverse_match,
                                                                     G.neighbors,
                                                                     H.neighbors,
                                                                     G.neighbors_ordered,
                                                                     H.neighbors_ordered)

        # evaluate semantic feasibility (always over original graphs)
        if(syntactic_feasibility):
            if(params.node_labels or params.edge_labels):
                semantic_feasibility = semantic_feasibility_undirected(params.node_labels,
                                                                       params.edge_labels,
                                                                       node1,
                                                                       node2,
                                                                       current_state.match_G,
                                                                       G.loops,
                                                                       current_state.forward_match,
                                                                       G.nodes,
                                                                       H.nodes,
                                                                       G.neighbors,
                                                                       G.edges,
                                                                       H.edges)

            # push to stack if valid
            if(semantic_feasibility):
                # extend match with the new pair
                change_in_state = extend_match_undirected(taken_from, node1, node2, G, H, current_state)

        # discard and finish if feasibility tests were not passed
        if((not syntactic_feasibility) or (not semantic_feasibility)):
            feasible_induced_anchor = False
            break

    # end of function
    return(feasible_induced_anchor)





# function: prepare initial state with induced anchor - directed ---------------
cdef cpp_bool search_isomorphisms_prepare_feasible_initial_state_directed(isomorphisms_search_params & params,
                                                                          isomorphisms_directed_graph & G,
                                                                          isomorphisms_directed_graph & H,
                                                                          isomorphisms_state_directed & current_state) noexcept:

    # local variables
    cdef cpp_bool inserted = False
    cdef cpp_bool semantic_feasibility = True
    cdef cpp_bool feasible_induced_anchor = True
    cdef cpp_bool in_syntactic_feasibility = True
    cdef cpp_bool out_syntactic_feasibility = True
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int taken_from = 0
    cdef cpp_pair[int, int] each_pair
    cdef isomorphisms_change_in_state_directed change_in_state

    # extend initial search state with each pair in the anchor
    for each_pair in params.induced_encoded_anchor:

        # split nodes
        node1 = each_pair.first
        node2 = each_pair.second

        # both vertices should be in the same ring or outside of them
        if(current_state.out_ring_G.find(node1) != current_state.out_ring_G.end()):
            taken_from = 0
            if(current_state.out_ring_H.find(node2) == current_state.out_ring_H.end()):
                feasible_induced_anchor = False
                break
        else:
            if(current_state.in_ring_G.find(node1) != current_state.in_ring_G.end()):
                taken_from = 1
                if(current_state.in_ring_H.find(node2) == current_state.in_ring_H.end()):
                    feasible_induced_anchor = False
                    break
            else:
                taken_from = 2
                if(current_state.out_ring_H.find(node2) != current_state.out_ring_H.end()):
                    feasible_induced_anchor = False
                    break
                if(current_state.in_ring_H.find(node2) != current_state.in_ring_H.end()):
                    feasible_induced_anchor = False
                    break

        # evaluate in syntactic feasibility (possibly over complement graph)
        if(params.complement):
            in_syntactic_feasibility = syntactic_feasibility_directed(node1,
                                                                      node2,
                                                                      current_state.in_ring_G,
                                                                      current_state.in_ring_H,
                                                                      current_state.out_ring_G,
                                                                      current_state.out_ring_H,
                                                                      current_state.match_G,
                                                                      current_state.match_H,
                                                                      G.loops,
                                                                      H.loops,
                                                                      current_state.forward_match,
                                                                      current_state.inverse_match,
                                                                      G.in_neighbors_complement,
                                                                      H.in_neighbors_complement,
                                                                      G.in_neighbors_ordered,
                                                                      H.in_neighbors_ordered)
        else:
            in_syntactic_feasibility = syntactic_feasibility_directed(node1,
                                                                      node2,
                                                                      current_state.in_ring_G,
                                                                      current_state.in_ring_H,
                                                                      current_state.out_ring_G,
                                                                      current_state.out_ring_H,
                                                                      current_state.match_G,
                                                                      current_state.match_H,
                                                                      G.loops,
                                                                      H.loops,
                                                                      current_state.forward_match,
                                                                      current_state.inverse_match,
                                                                      G.in_neighbors,
                                                                      H.in_neighbors,
                                                                      G.in_neighbors_ordered,
                                                                      H.in_neighbors_ordered)

        # evaluate out syntactic feasibility (possibly over complement graph)
        if(in_syntactic_feasibility):
            if(params.complement):
                out_syntactic_feasibility = syntactic_feasibility_directed(node1,
                                                                           node2,
                                                                           current_state.in_ring_G,
                                                                           current_state.in_ring_H,
                                                                           current_state.out_ring_G,
                                                                           current_state.out_ring_H,
                                                                           current_state.match_G,
                                                                           current_state.match_H,
                                                                           G.loops,
                                                                           H.loops,
                                                                           current_state.forward_match,
                                                                           current_state.inverse_match,
                                                                           G.out_neighbors_complement,
                                                                           H.out_neighbors_complement,
                                                                           G.out_neighbors_ordered,
                                                                           H.out_neighbors_ordered)
            else:
                out_syntactic_feasibility = syntactic_feasibility_directed(node1,
                                                                           node2,
                                                                           current_state.in_ring_G,
                                                                           current_state.in_ring_H,
                                                                           current_state.out_ring_G,
                                                                           current_state.out_ring_H,
                                                                           current_state.match_G,
                                                                           current_state.match_H,
                                                                           G.loops,
                                                                           H.loops,
                                                                           current_state.forward_match,
                                                                           current_state.inverse_match,
                                                                           G.out_neighbors,
                                                                           H.out_neighbors,
                                                                           G.out_neighbors_ordered,
                                                                           H.out_neighbors_ordered)

            # evaluate semantic feasibility (always over original graphs)
            if(out_syntactic_feasibility):
                if(params.node_labels or params.edge_labels):
                    semantic_feasibility = semantic_feasibility_directed(params.node_labels,
                                                                         params.edge_labels,
                                                                         node1,
                                                                         node2,
                                                                         current_state.match_G,
                                                                         G.loops,
                                                                         current_state.forward_match,
                                                                         G.nodes,
                                                                         H.nodes,
                                                                         G.in_neighbors,
                                                                         G.out_neighbors,
                                                                         G.edges,
                                                                         H.edges)

                # push to stack if valid
                if(semantic_feasibility):
                    # extend match with the new pair
                    change_in_state = extend_match_directed(taken_from, node1, node2, G, H, current_state)

        # discard and finish if feasibility tests were not passed
        if((not in_syntactic_feasibility) or (not out_syntactic_feasibility) or (not semantic_feasibility)):
            feasible_induced_anchor = False
            break

    # end of function
    return(feasible_induced_anchor)





# functions - search isomorphisms - undirected #################################





# function: core routine of VF2-like undirected approach - recursive -----------
cdef void search_isomorphisms_undirected(isomorphisms_search_params & params,
                                         isomorphisms_state_undirected & current_state,
                                         isomorphisms_undirected_graph & G,
                                         isomorphisms_undirected_graph & H,
                                         cpp_vector[cpp_set[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables
    cdef int node = 0
    cdef int matchable_node_G = 0
    cdef int candidate_node_H = 0
    cdef cpp_bool use_filter = False
    cdef cpp_bool viable_candidate = True
    cdef cpp_bool semantic_feasibility = True
    cdef cpp_bool syntactic_feasibility = True
    cdef cpp_pair[int, int] candidates_info
    cdef cpp_list[int] ordered_candidates
    cdef cpp_list[int] ring_G_ordered_backup
    cdef cpp_list[int] ring_H_ordered_backup
    cdef cpp_list[int] unmatched_G_ordered_backup
    cdef cpp_list[int] unmatched_H_ordered_backup
    cdef cpp_unordered_set[int] filtered_ring
    cdef isomorphisms_change_in_state_undirected change_in_state

    # save if isomorphism was reached
    if(current_state.match.size() == params.expected_order):
        all_matches.push_back(current_state.match)
        if(not params.all_isomorphisms):
            return

    # if not optimal yet then obtain available pairs
    if(current_state.match.size() < params.expected_order):

        # back-up state
        ring_G_ordered_backup = current_state.ring_G_ordered
        ring_H_ordered_backup = current_state.ring_H_ordered
        unmatched_G_ordered_backup = current_state.unmatched_G_ordered
        unmatched_H_ordered_backup = current_state.unmatched_H_ordered

        # get minimum and candidates to form pairs
        candidates_info = candidates_info_undirected(params, current_state)

        # get minimum node in G
        matchable_node_G = candidates_info.first

        # get list of candidates in H from ring
        if(candidates_info.second == 0):
            # order version of candidates
            ordered_candidates = current_state.ring_H_ordered
            # get filter for candidates, i.e., with compatibly matched neighbors in match (different from syntactic feasibility)
            filtered_ring = filter_candidates_undirected(matchable_node_G, G.neighbors_ordered, H.neighbors_ordered, current_state)
            # determine if filter is to be used
            use_filter = (filtered_ring.size() != ordered_candidates.size())
            # finish early if filter should be used but it is empty
            if(use_filter and (filtered_ring.empty())):
                return

        # get list of candidates in H from unmatched
        if(candidates_info.second == 1):
            # order version of candidates
            ordered_candidates = current_state.unmatched_H_ordered

        # evaluate candidate pairs
        for candidate_node_H in ordered_candidates:

            # test if candidate is in filter (only if required)
            if(use_filter):
                viable_candidate = (filtered_ring.find(candidate_node_H) != filtered_ring.end())

            # control variable changing only when filter is actually being used
            if(viable_candidate):

                # evaluate syntactic feasibility (possibly over complement graph)
                if(params.complement):
                    syntactic_feasibility = syntactic_feasibility_undirected(matchable_node_G,
                                                                             candidate_node_H,
                                                                             current_state.ring_G,
                                                                             current_state.ring_H,
                                                                             current_state.match_G,
                                                                             current_state.match_H,
                                                                             G.loops,
                                                                             H.loops,
                                                                             current_state.forward_match,
                                                                             current_state.inverse_match,
                                                                             G.neighbors_complement,
                                                                             H.neighbors_complement,
                                                                             G.neighbors_ordered,
                                                                             H.neighbors_ordered)
                else:
                    syntactic_feasibility = syntactic_feasibility_undirected(matchable_node_G,
                                                                             candidate_node_H,
                                                                             current_state.ring_G,
                                                                             current_state.ring_H,
                                                                             current_state.match_G,
                                                                             current_state.match_H,
                                                                             G.loops,
                                                                             H.loops,
                                                                             current_state.forward_match,
                                                                             current_state.inverse_match,
                                                                             G.neighbors,
                                                                             H.neighbors,
                                                                             G.neighbors_ordered,
                                                                             H.neighbors_ordered)

                # evaluate semantic feasibility (always over original graphs)
                if(syntactic_feasibility):
                    if(params.node_labels or params.edge_labels):
                        semantic_feasibility = semantic_feasibility_undirected(params.node_labels,
                                                                               params.edge_labels,
                                                                               matchable_node_G,
                                                                               candidate_node_H,
                                                                               current_state.match_G,
                                                                               G.loops,
                                                                               current_state.forward_match,
                                                                               G.nodes,
                                                                               H.nodes,
                                                                               G.neighbors,
                                                                               G.edges,
                                                                               H.edges)

                    # push to stack if valid
                    if(semantic_feasibility):
                        # extend match with the new pair
                        change_in_state = extend_match_undirected(candidates_info.second, matchable_node_G, candidate_node_H, G, H, current_state)
                        # recursive call
                        search_isomorphisms_undirected(params, current_state, G, H, all_matches)
                        # finish if only one isosmorphism was requested and it was already found
                        if(not all_matches.empty()):
                            if(not params.all_isomorphisms):
                                return
                        # restore state unordered sets
                        restore_match_undirected(matchable_node_G, candidate_node_H, change_in_state, current_state)
                        # restore state ordered lists
                        current_state.ring_G_ordered = ring_G_ordered_backup
                        current_state.ring_H_ordered = ring_H_ordered_backup
                        current_state.unmatched_G_ordered = unmatched_G_ordered_backup
                        current_state.unmatched_H_ordered = unmatched_H_ordered_backup

    # end of function





# function: get candidate pairs for undirected isomorphism search --------------
cdef cpp_pair[int, int] candidates_info_undirected(isomorphisms_search_params & params,
                                                   isomorphisms_state_undirected & current_state) noexcept:

    # output holders
    cdef cpp_pair[int, int] candidates_info

    # local variables
    cdef int node = 0
    cdef int minimum_node = -1
    cdef int minimum_value = 0

    # initialize with impossible minimum (total order ranges from 1 to expected_order)
    minimum_value = 10 + params.expected_order_int

    # NOTE: syntactic feasability for isomorphism guarantees that both rings have
    # always the same cardinality, thus we only have the following two disjoint options,
    # aided by the "filter" of the ring in H obtained in a second step.

    # build candidates either with nodes in the rings or with unmatched nodes
    if((not current_state.ring_G.empty()) and (not current_state.ring_H.empty())):

        # get node with minimum order from (ordered) ring in G
        for node in current_state.ring_G_ordered:
            if(params.total_order_G[node] < minimum_value):
                minimum_value = params.total_order_G[node]
                minimum_node = node

        # build output pair
        # NOTE: "filter" of candidates with compatibly-matched neighbors are obtained in a second step
        candidates_info.first = minimum_node
        candidates_info.second = 0

    else:

        # get node with minimum order from unmatched nodes in G
        minimum_node = current_state.unmatched_G_ordered.front()

        # build output pair
        candidates_info.first = minimum_node
        candidates_info.second = 1

    # end of function
    return(candidates_info)





# function: get compatibly-matched neighbors extending match -------------------
# NOTE: this is one of the "VF3" improvements
cdef cpp_unordered_set[int] filter_candidates_undirected(int matchable_node_G,
                                                         cpp_unordered_map[int, cpp_vector[int]] & neighbors_ordered_G,
                                                         cpp_unordered_map[int, cpp_vector[int]] & neighbors_ordered_H,
                                                         isomorphisms_state_undirected & current_state) noexcept:

    # output holders
    cdef cpp_unordered_set[int] filtered_candidates

    # local variables
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int mapped = 0

    # get filter
    for node1 in neighbors_ordered_G[matchable_node_G]:
        # continue as long as filter is not the whole ring
        if(filtered_candidates.size() != current_state.ring_H.size()):
            # only for neighbors in the match
            if(current_state.match_G.find(node1) != current_state.match_G.end()):
                # get image under match
                mapped = current_state.forward_match[node1]
                # traverse neighbors of image under match
                for node2 in neighbors_ordered_H[mapped]:
                    # save only the neighbors not in the match
                    if(current_state.match_H.find(node2) == current_state.match_H.end()):
                        # unordered set takes care of repetitions
                        filtered_candidates.insert(node2)
        else:
            # if the filter reached the whole ring then just finish early
            return(filtered_candidates)

    # end of function
    return(filtered_candidates)





# function: upgrade search state by extending match ----------------------------
cdef isomorphisms_change_in_state_undirected extend_match_undirected(int taken_from,
                                                                     int node1,
                                                                     int node2,
                                                                     isomorphisms_undirected_graph & G,
                                                                     isomorphisms_undirected_graph & H,
                                                                     isomorphisms_state_undirected & current_state) noexcept:

    # output holders
    cdef isomorphisms_change_in_state_undirected change_in_state

    # local variables
    cdef int node = 0
    cdef cpp_bool inserted
    cdef cpp_pair[int, int] temp_pair

    # reinitialize control variables
    change_in_state.removed_node_ring_G = False
    change_in_state.removed_node_ring_H = False
    inserted = False

    # add new pair
    temp_pair.first = node1
    temp_pair.second = node2
    current_state.match.insert(temp_pair)

    # add node to match in G
    current_state.match_G.insert(node1)

    # add node to match in H
    current_state.match_H.insert(node2)

    # upgrade forward map
    current_state.forward_match[node1] = node2

    # upgrade inverse map
    current_state.inverse_match[node2] = node1

    # remove node from unmatched nodes in G
    current_state.unmatched_G.erase(node1)
    current_state.unmatched_G_ordered.remove(node1)

    # remove node from unmatched nodes in H
    current_state.unmatched_H.erase(node2)
    current_state.unmatched_H_ordered.remove(node2)

    # remove node from ring in G if taken from there and the same for H
    if(taken_from == 0):
        # remove node from ring in G and assign true to change
        current_state.ring_G_ordered.remove(node1)
        change_in_state.removed_node_ring_G = current_state.ring_G.erase(node1)
        # remove node from ring in H and assign true to change
        current_state.ring_H_ordered.remove(node2)
        change_in_state.removed_node_ring_H = current_state.ring_H.erase(node2)

    # add unmatched neighbors to ring in G if not already there
    for node in G.neighbors_ordered[node1]:
        if(current_state.match_G.find(node) == current_state.match_G.end()):
            inserted = (current_state.ring_G.insert(node)).second
            if(inserted):
                current_state.ring_G_ordered.push_back(node)
                change_in_state.added_neighbors_ring_G.push(node)

    # add unmatched neighbors to ring in H if not already there
    for node in H.neighbors_ordered[node2]:
        if(current_state.match_H.find(node) == current_state.match_H.end()):
            inserted = (current_state.ring_H.insert(node)).second
            if(inserted):
                current_state.ring_H_ordered.push_back(node)
                change_in_state.added_neighbors_ring_H.push(node)

    # end of function
    return(change_in_state)





# function: backtrack search state by restoring match --------------------------
cdef void restore_match_undirected(int node1,
                                   int node2,
                                   isomorphisms_change_in_state_undirected & change_in_state,
                                   isomorphisms_state_undirected & current_state) noexcept:

    # local variables
    cdef int node = 0
    cdef cpp_pair[int, int] temp_pair

    # removed last added pair
    temp_pair.first = node1
    temp_pair.second = node2
    current_state.match.erase(temp_pair)

    # remove added node to match in G
    current_state.match_G.erase(node1)

    # remove added node to match in H
    current_state.match_H.erase(node2)

    # restore forward map
    current_state.forward_match.erase(node1)

    # restore inverse map
    current_state.inverse_match.erase(node2)

    # re-insert node to unnmatched nodes in G
    current_state.unmatched_G.insert(node1)

    # re-insert node to unnmatched nodes in H
    current_state.unmatched_H.insert(node2)

    # re-insert node to ring in G if taken from there
    if(change_in_state.removed_node_ring_G):
        current_state.ring_G.insert(node1)

    # re-insert node to ring in H if taken from there
    if(change_in_state.removed_node_ring_H):
        current_state.ring_H.insert(node2)

    # remove neighbors from ring in G if added to it
    while(not change_in_state.added_neighbors_ring_G.empty()):
        current_state.ring_G.erase(change_in_state.added_neighbors_ring_G.top())
        change_in_state.added_neighbors_ring_G.pop()

    # remove neighbors from ring in H if added to it
    while(not change_in_state.added_neighbors_ring_H.empty()):
        current_state.ring_H.erase(change_in_state.added_neighbors_ring_H.top())
        change_in_state.added_neighbors_ring_H.pop()

    # end of function





# function: evaluate syntactic feasability for isomorphism search --------------
cdef cpp_bool syntactic_feasibility_undirected(int node1,
                                               int node2,
                                               cpp_unordered_set[int] & ring_G,
                                               cpp_unordered_set[int] & ring_H,
                                               cpp_unordered_set[int] & current_match_G,
                                               cpp_unordered_set[int] & current_match_H,
                                               cpp_unordered_set[int] & loops_G,
                                               cpp_unordered_set[int] & loops_H,
                                               cpp_unordered_map[int, int] & forward_match,
                                               cpp_unordered_map[int, int] & inverse_match,
                                               cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                               cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                               cpp_unordered_map[int, cpp_vector[int]] & ordered_neigh_G,
                                               cpp_unordered_map[int, cpp_vector[int]] & ordered_neigh_H) noexcept:

    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef int neighbors_ring_G = 0
    cdef int neighbors_ring_H = 0
    cdef int neighbors_extern_G = 0
    cdef int neighbors_extern_H = 0

    # consistency of degree
    if(neigh_G[node1].size() != neigh_H[node2].size()):
        return(False)

    # loop-consistency-test
    if(loops_G.find(node1) != loops_G.end()):
        if(loops_H.find(node2) == loops_H.end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(loops_G.find(node1) == loops_G.end()):
        if(loops_H.find(node2) != loops_H.end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)

    # look ahead 0: consistency of neighbors in match, while doing tripartition of neighbors
    for node in ordered_neigh_G[node1]:
        if(current_match_G.find(node) != current_match_G.end()):
            # check that the mapping is also the corresponding neighbor
            mapped = forward_match[node]
            if(neigh_H[node2].find(mapped) == neigh_H[node2].end()):
                return(False)
        else:
            if(ring_G.find(node) != ring_G.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_G = neighbors_ring_G + 1
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_G = neighbors_extern_G + 1

    for node in ordered_neigh_H[node2]:
        if(current_match_H.find(node) != current_match_H.end()):
            # check that the mapping is also the corresponding neighbor
            mapped = inverse_match[node]
            if(neigh_G[node1].find(mapped) == neigh_G[node1].end()):
                return(False)
        else:
            if(ring_H.find(node) != ring_H.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_H = neighbors_ring_H + 1
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_H = neighbors_extern_H + 1

    # look ahead 1: consistency of neighbors in ring (not in match but adjacent to match)
    if(neighbors_ring_G != neighbors_ring_H):
        return(False)

    # look ahead 2: consistency of extern neighbors (neither in match nor adjacent to match)
    if(neighbors_extern_G != neighbors_extern_H):
        return(False)

    # end of function
    return(True)





# function: evaluate semantic feasability for isomorphism search ---------------
cdef cpp_bool semantic_feasibility_undirected(cpp_bool node_labels,
                                              cpp_bool edge_labels,
                                              int node1,
                                              int node2,
                                              cpp_unordered_set[int] & current_match_G,
                                              cpp_unordered_set[int] & loops_G,
                                              cpp_unordered_map[int, int] & forward_match,
                                              cpp_unordered_map[int, int] & nodes_G,
                                              cpp_unordered_map[int, int] & nodes_H,
                                              cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                              cpp_unordered_map[cpp_string, int] & edges_G,
                                              cpp_unordered_map[cpp_string, int] & edges_H) noexcept:

    # local variables
    cdef int node = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string labeled_edge_G
    cdef cpp_string labeled_edge_H

    # compare node-labels
    if(node_labels):
        if(nodes_G[node1] != nodes_H[node2]):
            return(False)

    # compare edge-labels
    if(edge_labels):

        # compare loop-labels
        if(loops_G.find(node1) != loops_G.end()):
            labeled_edge_G = to_string(node1) + comma + to_string(node1)
            labeled_edge_H = to_string(node2) + comma + to_string(node2)
            # compare labeled edges
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

        # compare non-loop labels
        for node in neigh_G[node1]:
            if(current_match_G.find(node) != current_match_G.end()):
                # edge in G with only one end in match
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                # edge in H with only one end in match
                labeled_edge_H = to_string(node2) + comma + to_string(forward_match[node])
                # compare labeled edges
                if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                    return(False)

    # end of function
    return(True)





# functions - search isomorphisms - directed ###################################





# function: core routine of VF2-like directed approach - recursive -------------
cdef void search_isomorphisms_directed(isomorphisms_search_params & params,
                                       isomorphisms_state_directed & current_state,
                                       isomorphisms_directed_graph & G,
                                       isomorphisms_directed_graph & H,
                                       cpp_vector[cpp_set[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables
    cdef int node = 0
    cdef int matchable_node_G = 0
    cdef int candidate_node_H = 0
    cdef cpp_bool use_filter = False
    cdef cpp_bool viable_candidate = True
    cdef cpp_bool semantic_feasibility = True
    cdef cpp_bool in_syntactic_feasibility = True
    cdef cpp_bool out_syntactic_feasibility = True
    cdef cpp_pair[int, int] candidates_info
    cdef cpp_list[int] ordered_candidates
    cdef cpp_list[int] in_ring_G_ordered_backup
    cdef cpp_list[int] out_ring_G_ordered_backup
    cdef cpp_list[int] in_ring_H_ordered_backup
    cdef cpp_list[int] out_ring_H_ordered_backup
    cdef cpp_list[int] unmatched_G_ordered_backup
    cdef cpp_list[int] unmatched_H_ordered_backup
    cdef cpp_unordered_set[int] filtered_ring
    cdef isomorphisms_change_in_state_directed change_in_state

    # save if isomorphism was reached
    if(current_state.match.size() == params.expected_order):
        all_matches.push_back(current_state.match)
        if(not params.all_isomorphisms):
            return

    # if not optimal yet then obtain available pairs
    if(current_state.match.size() < params.expected_order):

        # back-up state
        in_ring_G_ordered_backup = current_state.in_ring_G_ordered
        out_ring_G_ordered_backup = current_state.out_ring_G_ordered
        in_ring_H_ordered_backup = current_state.in_ring_H_ordered
        out_ring_H_ordered_backup = current_state.out_ring_H_ordered
        unmatched_G_ordered_backup = current_state.unmatched_G_ordered
        unmatched_H_ordered_backup = current_state.unmatched_H_ordered

        # get minimum and candidates to form pairs
        candidates_info = candidates_info_directed(params, current_state)

        # get minimum node in G
        matchable_node_G = candidates_info.first

        # get list of candidates in H from out-ring
        if(candidates_info.second == 0):
            # order version of candidates
            ordered_candidates = current_state.out_ring_H_ordered
            # get filter for candidates, i.e., with compatibly matched neighbors in match (different from syntactic feasibility)
            filtered_ring = filter_candidates_directed(matchable_node_G, G.in_neighbors_ordered, H.out_neighbors_ordered, current_state, current_state.out_ring_H)
            # determine if filter is to be used
            use_filter = (filtered_ring.size() != ordered_candidates.size())
            # finish early if filter should be used but it is empty
            if(use_filter and (filtered_ring.empty())):
                return

        # get list of candidates in H from in-ring
        if(candidates_info.second == 1):
            # order version of candidates
            ordered_candidates = current_state.in_ring_H_ordered
            # get filter for candidates, i.e., with compatibly matched neighbors in match (different from syntactic feasibility)
            filtered_ring = filter_candidates_directed(matchable_node_G, G.out_neighbors_ordered, H.in_neighbors_ordered, current_state, current_state.in_ring_H)
            # determine if filter is to be used
            use_filter = (filtered_ring.size() != ordered_candidates.size())
            # finish early if filter should be used but it is empty
            if(use_filter and (filtered_ring.empty())):
                return

        # get list of candidates in H from unmatched
        if(candidates_info.second == 2):
            # order version of candidates
            ordered_candidates = current_state.unmatched_H_ordered

        # evaluate candidate pairs
        for candidate_node_H in ordered_candidates:

            # test if candidate is in filter (only if required)
            if(use_filter):
                viable_candidate = (filtered_ring.find(candidate_node_H) != filtered_ring.end())

            # control variable changing only when filter is actually being used
            if(viable_candidate):

                # evaluate in syntactic feasibility (possibly over complement graph)
                if(params.complement):
                    in_syntactic_feasibility = syntactic_feasibility_directed(matchable_node_G,
                                                                              candidate_node_H,
                                                                              current_state.in_ring_G,
                                                                              current_state.in_ring_H,
                                                                              current_state.out_ring_G,
                                                                              current_state.out_ring_H,
                                                                              current_state.match_G,
                                                                              current_state.match_H,
                                                                              G.loops,
                                                                              H.loops,
                                                                              current_state.forward_match,
                                                                              current_state.inverse_match,
                                                                              G.in_neighbors_complement,
                                                                              H.in_neighbors_complement,
                                                                              G.in_neighbors_ordered,
                                                                              H.in_neighbors_ordered)
                else:
                    in_syntactic_feasibility = syntactic_feasibility_directed(matchable_node_G,
                                                                              candidate_node_H,
                                                                              current_state.in_ring_G,
                                                                              current_state.in_ring_H,
                                                                              current_state.out_ring_G,
                                                                              current_state.out_ring_H,
                                                                              current_state.match_G,
                                                                              current_state.match_H,
                                                                              G.loops,
                                                                              H.loops,
                                                                              current_state.forward_match,
                                                                              current_state.inverse_match,
                                                                              G.in_neighbors,
                                                                              H.in_neighbors,
                                                                              G.in_neighbors_ordered,
                                                                              H.in_neighbors_ordered)

                # evaluate out syntactic feasibility (possibly over complement graph)
                if(in_syntactic_feasibility):
                    if(params.complement):
                        out_syntactic_feasibility = syntactic_feasibility_directed(matchable_node_G,
                                                                                   candidate_node_H,
                                                                                   current_state.in_ring_G,
                                                                                   current_state.in_ring_H,
                                                                                   current_state.out_ring_G,
                                                                                   current_state.out_ring_H,
                                                                                   current_state.match_G,
                                                                                   current_state.match_H,
                                                                                   G.loops,
                                                                                   H.loops,
                                                                                   current_state.forward_match,
                                                                                   current_state.inverse_match,
                                                                                   G.out_neighbors_complement,
                                                                                   H.out_neighbors_complement,
                                                                                   G.out_neighbors_ordered,
                                                                                   H.out_neighbors_ordered)
                    else:
                        out_syntactic_feasibility = syntactic_feasibility_directed(matchable_node_G,
                                                                                   candidate_node_H,
                                                                                   current_state.in_ring_G,
                                                                                   current_state.in_ring_H,
                                                                                   current_state.out_ring_G,
                                                                                   current_state.out_ring_H,
                                                                                   current_state.match_G,
                                                                                   current_state.match_H,
                                                                                   G.loops,
                                                                                   H.loops,
                                                                                   current_state.forward_match,
                                                                                   current_state.inverse_match,
                                                                                   G.out_neighbors,
                                                                                   H.out_neighbors,
                                                                                   G.out_neighbors_ordered,
                                                                                   H.out_neighbors_ordered)

                    # evaluate semantic feasibility (always over original graphs)
                    if(out_syntactic_feasibility):
                        if(params.node_labels or params.edge_labels):
                            semantic_feasibility = semantic_feasibility_directed(params.node_labels,
                                                                                 params.edge_labels,
                                                                                 matchable_node_G,
                                                                                 candidate_node_H,
                                                                                 current_state.match_G,
                                                                                 G.loops,
                                                                                 current_state.forward_match,
                                                                                 G.nodes,
                                                                                 H.nodes,
                                                                                 G.in_neighbors,
                                                                                 G.out_neighbors,
                                                                                 G.edges,
                                                                                 H.edges)

                        # push to stack if valid
                        if(semantic_feasibility):
                            # extend match with the new pair
                            change_in_state = extend_match_directed(candidates_info.second, matchable_node_G, candidate_node_H, G, H, current_state)
                            # recursive call
                            search_isomorphisms_directed(params, current_state, G, H, all_matches)
                            # finish if only one isosmorphism was requested and it was already found
                            if(not all_matches.empty()):
                                if(not params.all_isomorphisms):
                                    return
                            # restore state unordered sets
                            restore_match_directed(matchable_node_G, candidate_node_H, change_in_state, current_state)
                            # restore state ordered lists
                            current_state.in_ring_G_ordered = in_ring_G_ordered_backup
                            current_state.out_ring_G_ordered = out_ring_G_ordered_backup
                            current_state.in_ring_H_ordered = in_ring_H_ordered_backup
                            current_state.out_ring_H_ordered = out_ring_H_ordered_backup
                            current_state.unmatched_G_ordered = unmatched_G_ordered_backup
                            current_state.unmatched_H_ordered = unmatched_H_ordered_backup

    # end of function





# function: get candidate pairs for directed isomorphism search ----------------
cdef cpp_pair[int, int] candidates_info_directed(isomorphisms_search_params & params,
                                                 isomorphisms_state_directed & current_state) noexcept:

    # output holders
    cdef cpp_pair[int, int] candidates_info

    # local variables
    cdef int node = 0
    cdef int minimum_node = -1
    cdef int minimum_value = 0

    # initialize with impossible minimum (total order ranges from 1 to expected_order)
    minimum_value = 10 + params.expected_order_int

    # NOTE: syntactic feasability for isomorphism guarantees that all the rings have
    # always the same cardinality, thus we only have the following three disjoint options,
    # aided by the "filter" of the ring in H obtained in a second step.

    # build candidates from out ring
    if((not current_state.out_ring_G.empty()) and (not current_state.out_ring_H.empty())):

        # get node with minimum order from out ring in G
        for node in current_state.out_ring_G_ordered:
            if(params.total_order_G[node] < minimum_value):
                minimum_value = params.total_order_G[node]
                minimum_node = node

        # build output pair
        # NOTE: "filter" of candidates with compatibly-matched neighbors are obtained in a second step
        candidates_info.first = minimum_node
        candidates_info.second = 0

    else:

        # build candidates from in ring
        if((not current_state.in_ring_G.empty()) and (not current_state.in_ring_H.empty())):

            # get node with minimum order from in ring in G
            for node in current_state.in_ring_G_ordered:
                if(params.total_order_G[node] < minimum_value):
                    minimum_value = params.total_order_G[node]
                    minimum_node = node

            # build output pair
            # NOTE: "filter" of candidates with compatibly-matched neighbors are obtained in a second step
            candidates_info.first = minimum_node
            candidates_info.second = 1

        else:

            # get node with minimum order from unmatched nodes in G
            minimum_node = current_state.unmatched_G_ordered.front()

            # build output pair
            candidates_info.first = minimum_node
            candidates_info.second = 2

    # end of function
    return(candidates_info)





# function: get compatibly-matched neighbors extending match -------------------
# NOTE: this is one of the "VF3" improvements
cdef cpp_unordered_set[int] filter_candidates_directed(int matchable_node_G,
                                                       cpp_unordered_map[int, cpp_vector[int]] & neighbors_ordered_G,
                                                       cpp_unordered_map[int, cpp_vector[int]] & neighbors_ordered_H,
                                                       isomorphisms_state_directed & current_state,
                                                       cpp_unordered_set[int] & ring_H) noexcept:

    # output holders
    cdef cpp_unordered_set[int] filtered_candidates

    # local variables
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int mapped = 0

    # get filter
    for node1 in neighbors_ordered_G[matchable_node_G]:
        # continue as long as filter is not the whole ring
        if(filtered_candidates.size() != ring_H.size()):
            # only for neighbors in the match
            if(current_state.match_G.find(node1) != current_state.match_G.end()):
                # get image under match
                mapped = current_state.forward_match[node1]
                # traverse neighbors of image under match
                for node2 in neighbors_ordered_H[mapped]:
                    # save only the neighbors not in the match
                    if(current_state.match_H.find(node2) == current_state.match_H.end()):
                        # unordered set takes care of repetitions
                        filtered_candidates.insert(node2)
        else:
            # if the filter reached the whole ring then just finish early
            return(filtered_candidates)

    # end of function
    return(filtered_candidates)





# function: upgrade search state by extending match ----------------------------
cdef isomorphisms_change_in_state_directed extend_match_directed(int taken_from,
                                                                 int node1,
                                                                 int node2,
                                                                 isomorphisms_directed_graph & G,
                                                                 isomorphisms_directed_graph & H,
                                                                 isomorphisms_state_directed & current_state) noexcept:

    # output holders
    cdef isomorphisms_change_in_state_directed change_in_state

    # local variables
    cdef int node = 0
    cdef cpp_bool inserted
    cdef cpp_pair[int, int] temp_pair

    # reinitialize control variables
    change_in_state.removed_node_in_ring_G = False
    change_in_state.removed_node_out_ring_G = False
    change_in_state.removed_node_in_ring_H = False
    change_in_state.removed_node_out_ring_H = False
    inserted = False

    # add new pair
    temp_pair.first = node1
    temp_pair.second = node2
    current_state.match.insert(temp_pair)

    # add node to match in G
    current_state.match_G.insert(node1)

    # add node to match in H
    current_state.match_H.insert(node2)

    # upgrade forward map
    current_state.forward_match[node1] = node2

    # upgrade inverse map
    current_state.inverse_match[node2] = node1

    # remove node from unmatched nodes in G
    current_state.unmatched_G.erase(node1)
    current_state.unmatched_G_ordered.remove(node1)

    # remove node from unmatched nodes in H
    current_state.unmatched_H.erase(node2)
    current_state.unmatched_H_ordered.remove(node2)

    # remove node from in-ring and/or out-ring in G if taken from there and the same for H
    if((taken_from == 0) or (taken_from == 1)):
        # try removing node from in-ring in G
        current_state.in_ring_G_ordered.remove(node1)
        change_in_state.removed_node_in_ring_G = current_state.in_ring_G.erase(node1)
        # try removing node from out-ring in G
        current_state.out_ring_G_ordered.remove(node1)
        change_in_state.removed_node_out_ring_G = current_state.out_ring_G.erase(node1)
        # try removing node from in-ring in H
        current_state.in_ring_H_ordered.remove(node2)
        change_in_state.removed_node_in_ring_H = current_state.in_ring_H.erase(node2)
        # try removing node from out-ring in H
        current_state.out_ring_H_ordered.remove(node2)
        change_in_state.removed_node_out_ring_H = current_state.out_ring_H.erase(node2)

    # add unmatched in-neighbors to in-ring in G if not already there
    for node in G.in_neighbors_ordered[node1]:
        if(current_state.match_G.find(node) == current_state.match_G.end()):
            inserted = (current_state.in_ring_G.insert(node)).second
            if(inserted):
                current_state.in_ring_G_ordered.push_back(node)
                change_in_state.added_neighbors_in_ring_G.push(node)

    # add unmatched out-neighbors to out-ring in G if not already there
    for node in G.out_neighbors_ordered[node1]:
        if(current_state.match_G.find(node) == current_state.match_G.end()):
            inserted = (current_state.out_ring_G.insert(node)).second
            if(inserted):
                current_state.out_ring_G_ordered.push_back(node)
                change_in_state.added_neighbors_out_ring_G.push(node)

    # add unmatched in-neighbors to in-ring in H if not already there
    for node in H.in_neighbors_ordered[node2]:
        if(current_state.match_H.find(node) == current_state.match_H.end()):
            inserted = (current_state.in_ring_H.insert(node)).second
            if(inserted):
                current_state.in_ring_H_ordered.push_back(node)
                change_in_state.added_neighbors_in_ring_H.push(node)

    # add unmatched out-neighbors to out-ring in H if not already there
    for node in H.out_neighbors_ordered[node2]:
        if(current_state.match_H.find(node) == current_state.match_H.end()):
            inserted = (current_state.out_ring_H.insert(node)).second
            if(inserted):
                current_state.out_ring_H_ordered.push_back(node)
                change_in_state.added_neighbors_out_ring_H.push(node)

    # end of function
    return(change_in_state)





# function: backtrack search state by restoring match --------------------------
cdef void restore_match_directed(int node1,
                                 int node2,
                                 isomorphisms_change_in_state_directed & change_in_state,
                                 isomorphisms_state_directed & current_state) noexcept:

    # local variables
    cdef int node = 0
    cdef cpp_pair[int, int] temp_pair

    # removed last added pair
    temp_pair.first = node1
    temp_pair.second = node2
    current_state.match.erase(temp_pair)

    # remove added node to match in G
    current_state.match_G.erase(node1)

    # remove added node to match in H
    current_state.match_H.erase(node2)

    # restore forward map
    current_state.forward_match.erase(node1)

    # restore inverse map
    current_state.inverse_match.erase(node2)

    # re-insert node to unnmatched nodes in G
    current_state.unmatched_G.insert(node1)

    # re-insert node to unnmatched nodes in H
    current_state.unmatched_H.insert(node2)

    # re-insert node to in-ring in G if taken from there
    if(change_in_state.removed_node_in_ring_G):
        current_state.in_ring_G.insert(node1)

    # re-insert node to out-ring in G if taken from there
    if(change_in_state.removed_node_out_ring_G):
        current_state.out_ring_G.insert(node1)

    # re-insert node to in-ring in H if taken from there
    if(change_in_state.removed_node_in_ring_H):
        current_state.in_ring_H.insert(node2)

    # re-insert node to out-ring in H if taken from there
    if(change_in_state.removed_node_out_ring_H):
        current_state.out_ring_H.insert(node2)

    # remove neighbors from in-ring in G if added to it
    while(not change_in_state.added_neighbors_in_ring_G.empty()):
        current_state.in_ring_G.erase(change_in_state.added_neighbors_in_ring_G.top())
        change_in_state.added_neighbors_in_ring_G.pop()

    # remove neighbors from out-ring in G if added to it
    while(not change_in_state.added_neighbors_out_ring_G.empty()):
        current_state.out_ring_G.erase(change_in_state.added_neighbors_out_ring_G.top())
        change_in_state.added_neighbors_out_ring_G.pop()

    # remove neighbors from in-ring in H if added to it
    while(not change_in_state.added_neighbors_in_ring_H.empty()):
        current_state.in_ring_H.erase(change_in_state.added_neighbors_in_ring_H.top())
        change_in_state.added_neighbors_in_ring_H.pop()

    # remove neighbors from out-ring in H if added to it
    while(not change_in_state.added_neighbors_out_ring_H.empty()):
        current_state.out_ring_H.erase(change_in_state.added_neighbors_out_ring_H.top())
        change_in_state.added_neighbors_out_ring_H.pop()

    # end of function





# function: evaluate syntactic feasability for isomorphism search --------------
cdef cpp_bool syntactic_feasibility_directed(int node1,
                                             int node2,
                                             cpp_unordered_set[int] & in_ring_G,
                                             cpp_unordered_set[int] & in_ring_H,
                                             cpp_unordered_set[int] & out_ring_G,
                                             cpp_unordered_set[int] & out_ring_H,
                                             cpp_unordered_set[int] & current_match_G,
                                             cpp_unordered_set[int] & current_match_H,
                                             cpp_unordered_set[int] & loops_G,
                                             cpp_unordered_set[int] & loops_H,
                                             cpp_unordered_map[int, int] & forward_match,
                                             cpp_unordered_map[int, int] & inverse_match,
                                             cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                             cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                             cpp_unordered_map[int, cpp_vector[int]] & ordered_neigh_G,
                                             cpp_unordered_map[int, cpp_vector[int]] & ordered_neigh_H) noexcept:

    # local variables
    cdef cpp_bool ring_neighbor = False
    cdef int node = 0
    cdef int mapped = 0
    cdef int neighbors_in_ring_G = 0
    cdef int neighbors_in_ring_H = 0
    cdef int neighbors_out_ring_G = 0
    cdef int neighbors_out_ring_H = 0
    cdef int neighbors_extern_G = 0
    cdef int neighbors_extern_H = 0

    # consistency of degrees
    if(neigh_G[node1].size() != neigh_H[node2].size()):
        return(False)

    # loop-consistency-test
    if(loops_G.find(node1) != loops_G.end()):
        if(loops_H.find(node2) == loops_H.end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(loops_G.find(node1) == loops_G.end()):
        if(loops_H.find(node2) != loops_H.end()):
            # node1 has no loop in G but node2 has a loop in H
            return(False)

    # look ahead 0: consistency of neighbors in match, while doing tripartition of neighbors
    for node in ordered_neigh_G[node1]:
        if(current_match_G.find(node) != current_match_G.end()):
            # check that the mapping is also the corresponding neighbor
            mapped = forward_match[node]
            if(neigh_H[node2].find(mapped) == neigh_H[node2].end()):
                return(False)
        else:
            # reinitialize flag
            ring_neighbor = False
            # neighbor in in-ring
            if(in_ring_G.find(node) != in_ring_G.end()):
                neighbors_in_ring_G = neighbors_in_ring_G + 1
                ring_neighbor = True
            # neighbor (possibly also) in out-ring
            if(out_ring_G.find(node) != out_ring_G.end()):
                neighbors_out_ring_G = neighbors_out_ring_G + 1
                ring_neighbor = True
            # neighbor not in any ring nor in match
            if(not ring_neighbor):
                neighbors_extern_G = neighbors_extern_G + 1

    for node in ordered_neigh_H[node2]:
        if(current_match_H.find(node) != current_match_H.end()):
            # check that the mapping is also the corresponding neighbor
            mapped = inverse_match[node]
            if(neigh_G[node1].find(mapped) == neigh_G[node1].end()):
                return(False)
        else:
            # reinitialize flag
            ring_neighbor = False
            # neighbor in in-ring
            if(in_ring_H.find(node) != in_ring_H.end()):
                neighbors_in_ring_H = neighbors_in_ring_H + 1
                ring_neighbor = True
            # neighbor (possibly also) in out-ring
            if(out_ring_H.find(node) != out_ring_H.end()):
                neighbors_out_ring_H = neighbors_out_ring_H + 1
                ring_neighbor = True
            # neighbor not in any ring nor in match
            if(not ring_neighbor):
                neighbors_extern_H = neighbors_extern_H + 1

    # look ahead 1: consistency of neighbors in rings (not in match but adjacent to match)
    if(neighbors_in_ring_G != neighbors_in_ring_H):
        return(False)
    if(neighbors_out_ring_G != neighbors_out_ring_H):
        return(False)

    # look ahead 2: consistency of extern neighbors (neither in match nor adjacent to match)
    if(neighbors_extern_G != neighbors_extern_H):
        return(False)

    # end of function
    return(True)





# function: evaluate semantic feasability for isomorphism search ---------------
cdef cpp_bool semantic_feasibility_directed(cpp_bool node_labels,
                                            cpp_bool edge_labels,
                                            int node1,
                                            int node2,
                                            cpp_unordered_set[int] & current_match_G,
                                            cpp_unordered_set[int] & loops_G,
                                            cpp_unordered_map[int, int] & forward_match,
                                            cpp_unordered_map[int, int] & nodes_G,
                                            cpp_unordered_map[int, int] & nodes_H,
                                            cpp_unordered_map[int, cpp_unordered_set[int]] & in_neigh_G,
                                            cpp_unordered_map[int, cpp_unordered_set[int]] & out_neigh_G,
                                            cpp_unordered_map[cpp_string, int] & edges_G,
                                            cpp_unordered_map[cpp_string, int] & edges_H) noexcept:

    # local variables
    cdef int node = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string labeled_edge_G
    cdef cpp_string labeled_edge_H

    # compare node-labels
    if(node_labels):
        if(nodes_G[node1] != nodes_H[node2]):
            return(False)

    # compare edge-labels
    if(edge_labels):

        # compare loop-labels
        if(loops_G.find(node1) != loops_G.end()):
            labeled_edge_G = to_string(node1) + comma + to_string(node1)
            labeled_edge_H = to_string(node2) + comma + to_string(node2)
            # compare labeled edges
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

        # compare non-loop edge labels with in-neighbors
        for node in in_neigh_G[node1]:
            if(current_match_G.find(node) != current_match_G.end()):
                # edge in G with in-neighbor in match
                labeled_edge_G = to_string(node) + comma + to_string(node1)
                # edge in H with out-neighbor in match
                labeled_edge_H = to_string(forward_match[node]) + comma + to_string(node2)
                # compare labeled edges
                if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                    return(False)

        # compare non-loop edge labels with out-neighbors
        for node in out_neigh_G[node1]:
            if(current_match_G.find(node) != current_match_G.end()):
                # edge in G with out-neighbor in match
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                # edge in H with out-neighbor in match
                labeled_edge_H = to_string(node2) + comma + to_string(forward_match[node])
                # compare labeled edges
                if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                    return(False)

    # end of function
    return(True)





################################################################################
################################################################################
