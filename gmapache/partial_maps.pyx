################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: partial_maps                                                       #
#                                                                              #
# - Description: analysis of properties of partial maps, like extensions into  #
#   isomorphisms, maximum induced anchored subgraphs, overlaps, and others.    #
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
from .integerization import encode_graphs, decode_graphs, encode_match, decode_match





# C/C++ structs ################################################################





# structs - graphs #############################################################





# struct: undirected graph -----------------------------------------------------
cdef struct partial_maps_undirected_graph:
    # node data
    cpp_unordered_map[int, int] nodes
    # degree sequence
    cpp_vector[int] degrees
    # loop data
    cpp_unordered_set[int] loops
    # edge data
    cpp_unordered_map[cpp_string, int] edges
    cpp_vector[cpp_pair[int, int]] raw_edges
    # neighbors data (for constant look-ups)
    cpp_unordered_map[int, cpp_unordered_set[int]] neighbors
    # ordered neighbors data (for insertion into ring)
    cpp_unordered_map[int, cpp_vector[int]] neighbors_ordered





# struct: directed graph -------------------------------------------------------
cdef struct partial_maps_directed_graph:
    # node data
    cpp_unordered_map[int, int] nodes
    # degree sequences
    cpp_vector[int] in_degrees
    cpp_vector[int] out_degrees
    # loop data
    cpp_unordered_set[int] loops
    # edge data
    cpp_unordered_map[cpp_string, int] edges
    cpp_vector[cpp_pair[int, int]] raw_edges
    # neighbors data (for constant look-ups)
    cpp_unordered_map[int, cpp_unordered_set[int]] in_neighbors
    cpp_unordered_map[int, cpp_unordered_set[int]] out_neighbors
    # ordered neighbors data (for insertion into ring)
    cpp_unordered_map[int, cpp_vector[int]] in_neighbors_ordered
    cpp_unordered_map[int, cpp_vector[int]] out_neighbors_ordered
    # neighbors in underlying undirected graph
    cpp_unordered_map[int, cpp_unordered_set[int]] connectivity_neighbors





# structs - states of search space  ############################################





# structure: state for the VF2-like partial maps search - undirected -----------
cdef struct partial_maps_state_undirected:
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





# structure: state for the VF2-like partial maps search - directed -------------
cdef struct partial_maps_state_directed:
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
cdef struct partial_maps_search_params:
    # flag from the caller wrapper
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs
    int caller
    # use node labels
    cpp_bool node_labels
    # use edge labels
    cpp_bool edge_labels
    # directed or undirected graphs
    cpp_bool directed_graphs
    # return all extensions
    cpp_bool all_extensions
    # search for matches with induced subgraph
    cpp_bool induced_subgraph
    # reachability for anchored-MCS search
    cpp_bool reachability
    # order of domain and codomain graphs
    int int_order_domain
    int int_order_codomain
    # expected amount of matches
    size_t expected_order
    int expected_order_int
    # total order for VF2-like search
    cpp_unordered_map[int, int] total_order_G
    cpp_unordered_map[int, int] total_order_H
    # inverse total order for VF2-like search
    cpp_unordered_map[int, int] inverse_total_order_G
    cpp_unordered_map[int, int] inverse_total_order_H
    # information of input anchor
    cpp_unordered_set[int] anchor_G
    cpp_unordered_set[int] anchor_H
    cpp_set[cpp_pair[int, int]] encoded_anchor
    # test for cover of degrees in unbalanced case
    cpp_bool test_unbalanced_degree_cover
    cpp_bool test_unbalanced_in_degree_cover
    cpp_bool test_unbalanced_out_degree_cover
    # backup data for parameters and other structures
    cpp_bool node_labels_backup
    cpp_bool edge_labels_backup
    cpp_unordered_set[int] unmatched_H_backup
    cpp_list[int] unmatched_H_ordered_backup





# structs - auxiliary structs for search #######################################





# struct: nodes added or removed from rings during recursive search ------------
cdef struct partial_maps_change_in_state_undirected:
    # node effectively removed from ring in G
    cpp_bool removed_node_ring_G
    # node effectively removed from ring in H
    cpp_bool removed_node_ring_H
    # neighbors effectively added to ring in G
    cpp_stack[int] added_neighbors_ring_G
    # neighbors effectively added to ring in H
    cpp_stack[int] added_neighbors_ring_H





# struct: nodes added or removed from rings during recursive search ------------
cdef struct partial_maps_change_in_state_directed:
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





# functions - search for complete induced extension - wrapper ##################





# function: callable wrapper for searching complete induced extension ----------
def search_complete_induced_extension(nx_G = nx.Graph(),           # can also be a networkx DiGraph
                                      nx_H = nx.Graph(),           # can also be a networkx DiGraph
                                      input_anchor = [],           # anchor partial map, should be a non-empty list
                                      node_labels = False,         # consider node labels when evaluating extension
                                      edge_labels = False,         # consider edge labels when evaluating extension
                                      all_extensions = False):     # by default stops when finding one extension (if any)

    # description
    """
    > description: receives two networkx graphs G and H of the same order, and a match
    between them (here called anchor), and uses a variant of the VF2 algorithm by doing
    an isomorphism-like search to obtain the induced complete extensions of the anchor,
    if any. These exist if and only if (1) the input graphs are "balanced" (with the
    same number of nodes and egdes of each label) and (2) if the provided partial map is
    a good atom map, i.e., if it already covers all the changing edges. Such extension,
    if it exists, is unique up to equivalence of atom maps, and by default the function
    will stop when finding it. This can be changed with the boolean parameter all_extensions,
    but it should be noted that such search can be more time consuming depending on the
    number of automotphisms of the input graphs. Such an extension, if it exists, is an
    isomorphism between the "reminder graphs", i.e., the graphs obtained from G and H
    when removing the reaction edges.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x
    from G and y from H. An exception is raised if the anchor is empty.
    * node_labels - boolean indicating if all node labels should be considered for the
    search or if they should be ignored (default).
    * edge_labels - boolean indicating if all edge labels should be considered for the
    search or if they should be ignored (default).
    * all_extensions - boolean indicating if the function should stop as soon as one
    complete extension is found (if any) -default behavior- or if it should search
    for all possible (complete) extensions. NOTE: mathematically speaking, the complete
    extension of a FIXED and GOOD anchor is unique up to equivalence of bijections,
    that is, up to equivalence of atom maps or correspondingly up to isomorphism of their
    ITS graphs. Nontheless it should be noted that changing the anchor may produce a
    non-equivalent (complete) extension. In other words, each call to this function can
    (mathematically) produce only one complete extension, but calls with different
    anchors can produce non-equivalent extensions, even if the anchors themselves produce
    isomorphic partial ITS graphs.

    > output:
    * extensions - (possibly empty) list of complete induced extensions, each as a list
    of 2-tuples (x, y) of nodes x from G and y from H representing the injective function
    preserving adjacency and (possibly) labels outside the anchor.
    * found_extensions - boolean value indicating if the extensions where found, i.e.,
    they cover all nodes of G and of H, and thus if they are bijections between G and H.
    If so, the anchor is what we refered to as a "good partial atom map", and equivalenteĺy
    the extension is an isomorphism between the reminder graphs obtained by removing the
    reaction edges from G and H.

    > calls:
    * gmapache.integerization.encode_graphs
    * gmapache.integerization.decode_graphs
    * gmapache.integerization.encode_match
    * gmapache.integerization.decode_match
    * gmapache.partial_maps.partial_maps_input_correctness
    * gmapache.partial_maps.partial_maps_order_nodes_concentric_reachability
    * gmapache.partial_maps.partial_maps_label_consistency
    * gmapache.partial_maps.partial_maps_order_neighbors
    * gmapache.partial_maps.partial_maps_prepare_initial_state_directed
    * gmapache.partial_maps.partial_maps_prepare_initial_state_undirected
    * gmapache.partial_maps.partial_maps_directed
    * gmapache.partial_maps.partial_maps_undirected
    """

    # output holders (python)
    found_extensions = False
    cdef list extensions = []

    # fixed threshold parameters
    cdef float scalation_value = 1.5

    # local variables (cython)
    cdef int i = 0
    cdef int deg = 0
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int counter = 0
    cdef int current_limit = 0
    cdef int required_limit = 0
    cdef cpp_bool reachable = True
    cdef cpp_bool input_correctness = True
    cdef cpp_bool consistent_labels = True
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_vector[int] deg_G
    cdef cpp_vector[int] deg_H
    cdef cpp_vector[int] in_deg_G
    cdef cpp_vector[int] in_deg_H
    cdef cpp_vector[int] out_deg_G
    cdef cpp_vector[int] out_deg_H
    cdef cpp_vector[int] temp_nodes
    cdef cpp_vector[int] each_vector
    cdef cpp_vector[int] all_nodes_G
    cdef cpp_vector[int] all_nodes_H
    cdef cpp_vector[cpp_vector[int]] all_edges_G
    cdef cpp_vector[cpp_vector[int]] all_edges_H
    cdef cpp_vector[cpp_set[cpp_pair[int, int]]] encoded_extensions
    cdef cpp_set[cpp_pair[int, int]] each_extension
    cdef cpp_unordered_map[int, int] node_degrees
    cdef partial_maps_search_params params
    cdef partial_maps_directed_graph directed_G
    cdef partial_maps_directed_graph directed_H
    cdef partial_maps_undirected_graph undirected_G
    cdef partial_maps_undirected_graph undirected_H
    cdef partial_maps_state_directed initial_state_directed
    cdef partial_maps_state_undirected initial_state_undirected

    # local variables (python)
    cdef list input_anchor_G = []
    cdef list input_anchor_H = []
    cdef list encoded_graphs = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_labels = dict()
    cdef dict encoded_edge_labels = dict()
    reachability = True
    x_obj = None
    y_obj = None
    node_obj = None
    nx_G_copy = None
    nx_H_copy = None
    undirected_copy = None

    # set caller parameter flag
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs
    params.caller = 0

    # test input correctness
    input_correctness = partial_maps_input_correctness(nx_G, nx_H, input_anchor, node_labels, edge_labels, all_extensions, reachability, params.caller)
    if(not input_correctness):
        return([], False)

    # quick test by comparing the degree sequences after removing the nodes in the anchor
    params.directed_graphs = nx.is_directed(nx_G)
    # copy input graphs
    nx_G_copy = deepcopy(nx_G)
    nx_H_copy = deepcopy(nx_H)
    # remove vertices in the anchor
    input_anchor_G = [x_obj for (x_obj, y_obj) in input_anchor]
    input_anchor_H = [y_obj for (x_obj, y_obj) in input_anchor]
    nx_G_copy.remove_nodes_from(input_anchor_G)
    nx_H_copy.remove_nodes_from(input_anchor_H)
    # obtain and compare degree sequences
    if(params.directed_graphs):
        # get in-degrees
        in_deg_G = [deg for (node_obj, deg) in list(nx_G_copy.in_degree())]
        in_deg_H = [deg for (node_obj, deg) in list(nx_H_copy.in_degree())]
        # sort in-degrees
        sort(in_deg_G.begin(), in_deg_G.end())
        sort(in_deg_H.begin(), in_deg_H.end())
        # compare in-degree sequences
        if(not in_deg_G == in_deg_H):
            return([], False)
        # get out-degrees
        out_deg_G = [deg for (node_obj, deg) in list(nx_G_copy.out_degree())]
        out_deg_H = [deg for (node_obj, deg) in list(nx_H_copy.out_degree())]
        # sort out-degrees
        sort(out_deg_G.begin(), out_deg_G.end())
        sort(out_deg_H.begin(), out_deg_H.end())
        # compare out-degree sequences
        if(not out_deg_G == out_deg_H):
            return([], False)
    else:
        # get degrees
        deg_G = [deg for (node_obj, deg) in list(nx_G_copy.degree())]
        deg_H = [deg for (node_obj, deg) in list(nx_H_copy.degree())]
        # sort degrees
        sort(deg_G.begin(), deg_G.end())
        sort(deg_H.begin(), deg_H.end())
        # compare degree sequences
        if(not deg_G == deg_H):
            return([], False)

    # save input parameters
    params.node_labels = node_labels
    params.edge_labels = edge_labels
    params.all_extensions = all_extensions
    params.expected_order = nx_G.order()
    params.expected_order_int = nx_G.order()
    params.reachability = True
    params.induced_subgraph = True

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_labels, encoded_edge_labels = encode_graphs([nx_G, nx_H])

    # encode match
    params.encoded_anchor = encode_match(input_anchor, encoded_node_names)
    for each_pair in params.encoded_anchor:
        params.anchor_G.insert(each_pair.first)
        params.anchor_H.insert(each_pair.second)

    # prepare and order nodes in concentric sets around the anchor
    if(params.directed_graphs):
        # prepare information for ordering of G
        temp_nodes = list(encoded_graphs[0].nodes())
        undirected_copy = deepcopy(encoded_graphs[0])
        undirected_copy = undirected_copy.to_undirected()
        directed_G.connectivity_neighbors = {node:set(undirected_copy.neighbors(node)) for node in temp_nodes}
        directed_G.nodes = {node:encoded_graphs[0].nodes[node]["GMNL"] for node in temp_nodes}
        # out degrees for ordering
        node_degrees = {node:encoded_graphs[0].out_degree(node) for node in temp_nodes}
        # order for directed G
        reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_G,
                                                                     temp_nodes,
                                                                     directed_G.connectivity_neighbors,
                                                                     node_degrees,
                                                                     all_nodes_G)
        # prepare information for ordering of H
        temp_nodes = list(encoded_graphs[1].nodes())
        undirected_copy = deepcopy(encoded_graphs[1])
        undirected_copy = undirected_copy.to_undirected()
        directed_H.connectivity_neighbors = {node:set(undirected_copy.neighbors(node)) for node in temp_nodes}
        directed_H.nodes = {node:encoded_graphs[1].nodes[node]["GMNL"] for node in temp_nodes}
        # out degrees for ordering
        node_degrees = {node:encoded_graphs[1].out_degree(node) for node in temp_nodes}
        # order for directed H
        reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_H,
                                                                     temp_nodes,
                                                                     directed_H.connectivity_neighbors,
                                                                     node_degrees,
                                                                     all_nodes_H)
    else:
        # prepare information for ordering of G
        temp_nodes = list(encoded_graphs[0].nodes())
        undirected_G.neighbors = {node:set(encoded_graphs[0].neighbors(node)) for node in temp_nodes}
        undirected_G.nodes = {node:encoded_graphs[0].nodes[node]["GMNL"] for node in temp_nodes}
        # degrees for ordering
        node_degrees = {node:encoded_graphs[0].degree(node) for node in temp_nodes}
        # order for undirected G
        reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_G,
                                                                     temp_nodes,
                                                                     undirected_G.neighbors,
                                                                     node_degrees,
                                                                     all_nodes_G)
        # prepare information for ordering of H
        temp_nodes = list(encoded_graphs[1].nodes())
        undirected_H.neighbors = {node:set(encoded_graphs[1].neighbors(node)) for node in temp_nodes}
        undirected_H.nodes = {node:encoded_graphs[1].nodes[node]["GMNL"] for node in temp_nodes}
        # degrees for ordering
        node_degrees = {node:encoded_graphs[1].degree(node) for node in temp_nodes}
        # order for undirected H
        reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_H,
                                                                     temp_nodes,
                                                                     undirected_H.neighbors,
                                                                     node_degrees,
                                                                     all_nodes_H)

    # prepare total orders
    if(params.directed_graphs):
        # nodes of directed G
        directed_G.loops = list(nx.nodes_with_selfloops(encoded_graphs[0]))
        counter = 0
        for node in all_nodes_G:
            # neighbors for G
            directed_G.in_neighbors[node] = set(encoded_graphs[0].predecessors(node))
            directed_G.out_neighbors[node] = set(encoded_graphs[0].neighbors(node))
            # save total order for the node
            counter = counter + 1
            params.total_order_G[node] = counter
            params.inverse_total_order_G[counter] = node
            initial_state_directed.unmatched_G.insert(node)
            initial_state_directed.unmatched_G_ordered.push_back(node)
        # nodes of directed H
        directed_H.loops = list(nx.nodes_with_selfloops(encoded_graphs[1]))
        counter = 0
        for node in all_nodes_H:
            # neighbors for H
            directed_H.in_neighbors[node] = set(encoded_graphs[1].predecessors(node))
            directed_H.out_neighbors[node] = set(encoded_graphs[1].neighbors(node))
            # save total order for the node
            counter = counter + 1
            params.total_order_H[node] = counter
            params.inverse_total_order_H[counter] = node
            initial_state_directed.unmatched_H.insert(node)
            initial_state_directed.unmatched_H_ordered.push_back(node)
    else:
        # nodes of undirected G
        undirected_G.loops = list(nx.nodes_with_selfloops(encoded_graphs[0]))
        counter = 0
        for node in all_nodes_G:
            # save total order for the node
            counter = counter + 1
            params.total_order_G[node] = counter
            params.inverse_total_order_G[counter] = node
            initial_state_undirected.unmatched_G.insert(node)
            initial_state_undirected.unmatched_G_ordered.push_back(node)
        # nodes of undirected H
        undirected_H.loops = list(nx.nodes_with_selfloops(encoded_graphs[1]))
        counter = 0
        for node in all_nodes_H:
            # save total order for the node
            counter = counter + 1
            params.total_order_H[node] = counter
            params.inverse_total_order_H[counter] = node
            initial_state_undirected.unmatched_H.insert(node)
            initial_state_undirected.unmatched_H_ordered.push_back(node)

    # prepare edges
    if(params.edge_labels):
        if(params.directed_graphs):
            # edges for G
            all_edges_G = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[0].edges(data = True)]
            for each_vector in all_edges_G:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                directed_G.raw_edges.push_back(temp_pair)
                # save labeled edge
                temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                directed_G.edges[temp_str] = each_vector[2]
            # edges for H
            all_edges_H = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[1].edges(data = True)]
            for each_vector in all_edges_H:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                directed_H.raw_edges.push_back(temp_pair)
                # save labeled edge
                temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                directed_H.edges[temp_str] = each_vector[2]
        else:
            # edges for G
            all_edges_G = [(node1, node2, info["GMEL"]) for (node1, node2, info) in encoded_graphs[0].edges(data = True)]
            for each_vector in all_edges_G:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                undirected_G.raw_edges.push_back(temp_pair)
                # save labeled edge
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
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                undirected_H.raw_edges.push_back(temp_pair)
                # save labeled edge
                if(each_vector[0] == each_vector[1]):
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[0])
                    undirected_H.edges[temp_str] = each_vector[2]
                else:
                    # save the two label edges to simplify future access
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                    undirected_H.edges[temp_str] = each_vector[2]
                    temp_str = to_string(each_vector[1]) + comma + to_string(each_vector[0])
                    undirected_H.edges[temp_str] = each_vector[2]

    # test consistency of labels if required
    if(params.node_labels or params.edge_labels):
        # test consistency of node and/or edge labels
        if(params.directed_graphs):
            consistent_labels = partial_maps_label_consistency(params.caller, params.node_labels, params.edge_labels,
                                                               directed_G.nodes, directed_H.nodes,
                                                               directed_G.edges, directed_H.edges,
                                                               directed_G.raw_edges, directed_H.raw_edges,
                                                               params.anchor_G, params.anchor_H)
        else:
            consistent_labels = partial_maps_label_consistency(params.caller, params.node_labels, params.edge_labels,
                                                               undirected_G.nodes, undirected_H.nodes,
                                                               undirected_G.edges, undirected_H.edges,
                                                               undirected_G.raw_edges, undirected_H.raw_edges,
                                                               params.anchor_G, params.anchor_H)
        # if labels are not consistent then the reminder graphs cannot be isomorphic
        if(not consistent_labels):
            return([], False)

    # prepare ordered neighbors
    if(params.directed_graphs):
        # order nodes of G
        partial_maps_order_neighbors(all_nodes_G, directed_G.in_neighbors, directed_G.in_neighbors_ordered,
                                     params.total_order_G, params.inverse_total_order_G)
        partial_maps_order_neighbors(all_nodes_G, directed_G.out_neighbors, directed_G.out_neighbors_ordered,
                                     params.total_order_G, params.inverse_total_order_G)
        # order nodes of H
        partial_maps_order_neighbors(all_nodes_H, directed_H.in_neighbors, directed_H.in_neighbors_ordered,
                                     params.total_order_H, params.inverse_total_order_H)
        partial_maps_order_neighbors(all_nodes_H, directed_H.out_neighbors, directed_H.out_neighbors_ordered,
                                     params.total_order_H, params.inverse_total_order_H)
    else:
        # order nodes of G
        partial_maps_order_neighbors(all_nodes_G, undirected_G.neighbors, undirected_G.neighbors_ordered,
                                     params.total_order_G, params.inverse_total_order_G)
        # order nodes of H
        partial_maps_order_neighbors(all_nodes_H, undirected_H.neighbors, undirected_H.neighbors_ordered,
                                     params.total_order_H, params.inverse_total_order_H)

    # set recursion limit
    required_limit = nx_G.order()
    current_limit = getrecursionlimit()
    if(current_limit < (scalation_value * required_limit)):
        setrecursionlimit(int(scalation_value * required_limit))

    # prepare initial state updated with anchor
    if(params.directed_graphs):
        partial_maps_prepare_initial_state_directed(params.encoded_anchor, params.anchor_G, params.anchor_H, directed_G, directed_H, initial_state_directed)
    else:
        partial_maps_prepare_initial_state_undirected(params.encoded_anchor, params.anchor_G, params.anchor_H, undirected_G, undirected_H, initial_state_undirected)

    # evaluate complete induced extension
    if(params.directed_graphs):
        partial_maps_directed(params, initial_state_directed, directed_G, directed_H, encoded_extensions)
    else:
        partial_maps_undirected(params, initial_state_undirected, undirected_G, undirected_H, encoded_extensions)

    # decode extensions
    for each_extension in encoded_extensions:
        extensions.append(decode_match(list(each_extension), encoded_node_names))

    # check if there were any extensions
    if(len(extensions) > 0):
        found_extensions = True

    # end of function
    return(extensions, found_extensions)





# functions - search maximum common anchored subgraphs - wrapper ###############





# function: callable wrapper for maximum common induced anchored subgraphs -----
def search_maximum_common_anchored_subgraphs(nx_G = nx.Graph(),           # can also be a networkx DiGraph
                                             nx_H = nx.Graph(),           # can also be a networkx DiGraph
                                             input_anchor = [],           # anchor partial map, should be a non-empty list
                                             node_labels = False,         # consider node labels when evaluating extension
                                             edge_labels = False,         # consider edge labels when evaluating extension
                                             all_extensions = False,      # by default stops when finding one extension (if any)
                                             reachability = True):        # all nodes should be reachable from at least one anchor node

    # description
    """
    > description: receives two non-null networkx (di-) graphs G and H (possibly with different number
    of nodes), and a non-empty injective map between them (here called anchor), and uses a variant
    of the VF2 algorithm by doing an isomorphism-like search and an iterative trimming to obtain
    the maximum common induced subgraphs that extend the matches in the anchor. The function
    specifically searches for proper extensions of the anchor, i.e., containing at least one more
    match than the anchor itself. If no proper extension is found then the function returns the anchor
    as a list and a boolean variable with value False. The parameter "reachability", controls if the
    candidate common subgraphs should be "connected" to the anchor, or better put "reachable" from
    the anchor in the sense that such subgraphs should contain a path between every node and at least
    one node from the anchor. Thus, if the anchor is connected and reachability is set to True, the
    function will equivalently get a maximum-common-induced-connected-anchored-subgraph. If both
    graphs have the same order, the function first will search for a complete-induced-extension and
    only if no such extension is found it will continue with the search for the maximum common induced
    anchored subgraphs, thus this function can be more time consuming that simply running the search
    for the complete induced extension. Moreover, if both graphs have the same order and a complete
    induced extension exists between them, this function will return such extension independently of
    the parameter "reachability" and regardless if the complete extension induces a connected ITS.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * input_anchor - inyective map as a non-empty list of 2-tuples (x, y) of nodes x from G
    and y from H. An exception is raised if the anchor is empty.
    * node_labels - boolean indicating if all labels of nodes (outside the anchor) should be
    considered for the search or if all of them should be ignored (default).
    * edge_labels - boolean indicating if all labels of edges (with at least one end outside the
    anchor) should be considered for the search or if all of them should be ignored (default).
    * all_extensions - boolean indicating if the function should stop as soon as one maximum
    induced subgraph is found (if any) -default behavior- or if it should search for all such
    possible common subgraphs extending the anchor. Here maximum refers to number of nodes.
    * reachability - boolean varibale that controls whether the nodes in the candidate subgraph
    should contain at least one path between each node and at least one anchor node (default).
    If the anchor is connected and reachability is set to True, then the search will return a
    maximum-common-induced-connected-subgraph properly containing the anchor, if any.

    > output:
    * extensions - non-empty list of maximum common induced subgraphs extending the anchor, each
    as a list of 2-tuples (x, y) of nodes x from G and y from H representing the injective function
    preserving adjacency and also labels outside the anchor, if required. If no proper extension
    was found then the anchor itself is returned inside this list.
    * found_proper_extensions - boolean value indicating if any proper maximum extensions of the
    anchor were found, i.e., containing the matches from the anchor but strictly bigger than it.

    > calls:
    * gmapache.integerization.encode_graphs
    * gmapache.integerization.decode_graphs
    * gmapache.integerization.encode_match
    * gmapache.integerization.decode_match
    * gmapache.partial_maps.partial_maps_input_correctness
    * gmapache.partial_maps.search_complete_induced_extension
    * gmapache.partial_maps.partial_maps_order_nodes_concentric_reachability
    * gmapache.partial_maps.partial_maps_order_neighbors
    * gmapache.partial_maps.partial_maps_iterative_trimming
    """

    # output holders (python)
    found_proper_extensions = False
    cdef list extensions = []

    # fixed threshold parameters
    cdef float scalation_value = 1.5

    # local variables (cython)
    cdef int i = 0
    cdef int deg = 0
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef int offset = 0
    cdef int covered = 0
    cdef int counter = 0
    cdef int order_G = 0
    cdef int order_H = 0
    cdef int order_bigger = 0
    cdef int order_smaller = 0
    cdef int current_limit = 0
    cdef int required_limit = 0
    cdef cpp_bool twisted = False
    cdef cpp_bool reachable = True
    cdef cpp_bool input_correctness = True
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_vector[int] deg_bigger
    cdef cpp_vector[int] deg_smaller
    cdef cpp_vector[int] in_deg_bigger
    cdef cpp_vector[int] in_deg_smaller
    cdef cpp_vector[int] out_deg_bigger
    cdef cpp_vector[int] out_deg_smaller
    cdef cpp_vector[int] temp_nodes
    cdef cpp_vector[int] all_nodes_bigger
    cdef cpp_vector[cpp_vector[int]] all_edges_bigger
    cdef cpp_vector[cpp_vector[int]] all_edges_smaller
    cdef cpp_vector[cpp_set[cpp_pair[int, int]]] encoded_extensions
    cdef cpp_vector[cpp_set[cpp_pair[int, int]]] untwisted_extensions
    cdef cpp_set[cpp_pair[int, int]] each_extension
    cdef cpp_set[cpp_pair[int, int]] temp_extension
    cdef cpp_set[cpp_pair[int, int]] original_anchor_encoded
    cdef cpp_unordered_map[int, int] node_degrees
    cdef cpp_unordered_map[int, int] all_removable_nodes
    cdef partial_maps_search_params params
    cdef partial_maps_directed_graph directed_bigger
    cdef partial_maps_directed_graph directed_smaller
    cdef partial_maps_undirected_graph undirected_bigger
    cdef partial_maps_undirected_graph undirected_smaller

    # local variables (python)
    cdef list encoded_graphs = []
    cdef dict info = dict()
    cdef dict encoded_node_names = dict()
    cdef dict encoded_node_labels = dict()
    cdef dict encoded_edge_labels = dict()
    node_obj = None
    nx_bigger = None
    nx_smaller = None
    undirected_copy = None

    # set caller parameter flag
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs
    params.caller = 1

    # test input correctness
    input_correctness = partial_maps_input_correctness(nx_G, nx_H, input_anchor, node_labels, edge_labels, all_extensions, reachability, params.caller)
    if(not input_correctness):
        return([input_anchor], False)

    # get order of graphs for local processing
    order_G = nx_G.order()
    order_H = nx_H.order()

    # only if the graphs have the same order we test for complete extension
    # NOTE: search_complete_induced_extension evaluates consistency of degree sequences
    if(order_G == order_H):
        extensions, found_proper_extensions = search_complete_induced_extension(nx_G = nx_G,
                                                                                nx_H = nx_H,
                                                                                input_anchor = input_anchor,
                                                                                node_labels = node_labels,
                                                                                edge_labels = edge_labels,
                                                                                all_extensions = all_extensions)
        # only if extensions were found then return, otherwise we continue with the normal search
        if(found_proper_extensions):
            return(extensions, found_proper_extensions)

    # save input parameters
    params.directed_graphs = nx.is_directed(nx_G)
    params.node_labels_backup = node_labels
    params.edge_labels_backup = edge_labels
    params.all_extensions = all_extensions
    params.reachability = reachability
    params.induced_subgraph = True
    params.test_unbalanced_degree_cover = False
    params.test_unbalanced_in_degree_cover = False
    params.test_unbalanced_out_degree_cover = False

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_labels, encoded_edge_labels = encode_graphs([nx_G, nx_H])

    # determine aliases for the smaller and bigger graph
    if(params.directed_graphs):
        if(order_G <= order_H):
            # set twisted flag
            twisted = False
            # take graphs preserving input order
            nx_smaller = encoded_graphs[0]
            nx_bigger = encoded_graphs[1]
        else:
            # set twisted flag
            twisted = True
            # take graphs inverting input order
            nx_smaller = encoded_graphs[1]
            nx_bigger = encoded_graphs[0]
    else:
        if(order_G <= order_H):
            # set twisted flag
            twisted = False
            # take graphs preserving input order
            nx_smaller = encoded_graphs[0]
            nx_bigger = encoded_graphs[1]
        else:
            # set twisted flag
            twisted = True
            # take graphs inverting input order
            nx_smaller = encoded_graphs[1]
            nx_bigger = encoded_graphs[0]

    # get order of graphs under aliases
    if(twisted):
        order_smaller = nx_H.order()
        order_bigger = nx_G.order()
        params.int_order_codomain = nx_G.order()
    else:
        order_smaller = nx_G.order()
        order_bigger = nx_H.order()
        params.int_order_codomain = nx_H.order()

    # test for cover of degree sequences between original graphs. If the cover doesnt hold then
    # we test for the cover in each candidate subgraph produced by the iterative trimming
    if(params.directed_graphs):
        # get in-degrees
        in_deg_smaller = [deg for (node_obj, deg) in list(nx_smaller.in_degree())]
        in_deg_bigger = [deg for (node_obj, deg) in list(nx_bigger.in_degree())]
        # sort in-degrees
        sort(in_deg_smaller.begin(), in_deg_smaller.end())
        sort(in_deg_bigger.begin(), in_deg_bigger.end())
        # save in-degree sequence of bigger for future references
        directed_smaller.in_degrees = in_deg_smaller
        directed_bigger.in_degrees = in_deg_bigger
        # test if in-degrees of bigger cover in-degrees of smaller
        i = 0
        offset = 0
        covered = 0
        for deg in in_deg_smaller:
            # decide if continue searching
            if(offset >= order_bigger):
                break
            else:
                # search for covering degree
                for i in range(offset, order_bigger):
                    if(deg <= in_deg_bigger[i]):
                        covered = covered + 1
                        break
                # get next offset
                offset = i + 1
        # test for in-cover
        if(covered < order_smaller):
            params.test_unbalanced_in_degree_cover = True
        # get out-degrees
        out_deg_smaller = [deg for (node_obj, deg) in list(nx_smaller.out_degree())]
        out_deg_bigger = [deg for (node_obj, deg) in list(nx_bigger.out_degree())]
        # sort out-degrees
        sort(out_deg_smaller.begin(), out_deg_smaller.end())
        sort(out_deg_bigger.begin(), out_deg_bigger.end())
        # save out-degree sequence of bigger for future references
        directed_smaller.out_degrees = out_deg_smaller
        directed_bigger.out_degrees = out_deg_bigger
        # test if out-degrees of bigger cover out-degrees of smaller
        i = 0
        offset = 0
        covered = 0
        for deg in out_deg_smaller:
            # decide if continue searching
            if(offset >= order_bigger):
                break
            else:
                # search for covering degree
                for i in range(offset, order_bigger):
                    if(deg <= out_deg_bigger[i]):
                        covered = covered + 1
                        break
                # get next offset
                offset = i + 1
        # test for out-cover
        if(covered < order_smaller):
            params.test_unbalanced_out_degree_cover = True
    else:
        # get degrees
        deg_smaller = [deg for (node_obj, deg) in list(nx_smaller.degree())]
        deg_bigger = [deg for (node_obj, deg) in list(nx_bigger.degree())]
        # sort degrees
        sort(deg_smaller.begin(), deg_smaller.end())
        sort(deg_bigger.begin(), deg_bigger.end())
        # save degree sequence of bigger for future references
        undirected_smaller.degrees = deg_smaller
        undirected_bigger.degrees = deg_bigger
        # test if degrees of bigger cover degrees of smaller
        i = 0
        offset = 0
        covered = 0
        for deg in deg_smaller:
            # decide if continue searching
            if(offset >= order_bigger):
                break
            else:
                # search for covering degree
                for i in range(offset, order_bigger):
                    if(deg <= deg_bigger[i]):
                        covered = covered + 1
                        break
                # get next offset
                offset = i + 1
        # test for cover of all degrees
        if(covered < order_smaller):
            params.test_unbalanced_degree_cover = True

    # encode anchor preserving original order as given in the input
    original_anchor_encoded = encode_match(input_anchor, encoded_node_names)

    # save encoded anchor into search parameters inverting input order if necessary
    # NOTE: inside the intensive routines and in associated search-parameters, the
    # domain graph is called G, while the codomain graph is always called H.
    if(twisted):
        # save inverted anchor
        for each_pair in original_anchor_encoded:
            temp_pair.first = each_pair.second
            temp_pair.second = each_pair.first
            params.anchor_G.insert(temp_pair.first)
            params.anchor_H.insert(temp_pair.second)
            params.encoded_anchor.insert(temp_pair)
    else:
        # save anchor
        params.encoded_anchor = original_anchor_encoded
        # save anchor in each graph
        for each_pair in params.encoded_anchor:
            params.anchor_G.insert(each_pair.first)
            params.anchor_H.insert(each_pair.second)

    # prepare nodes and neighbors for both graphs and ordered nodes only for bigger graph
    # NOTE: inside the intensive routines and in associated search-parameters, the
    # domain graph is called G, while the codomain graph is always called H.
    if(params.directed_graphs):
        # prepare nodes and neighbors of smaller directed graph
        directed_smaller.loops = list(nx.nodes_with_selfloops(nx_smaller))
        temp_nodes = list(nx_smaller.nodes())
        directed_smaller.nodes = {node:nx_smaller.nodes[node]["GMNL"] for node in temp_nodes}
        undirected_copy = deepcopy(nx_smaller)
        undirected_copy = undirected_copy.to_undirected()
        directed_smaller.connectivity_neighbors = {node:set(undirected_copy.neighbors(node)) for node in temp_nodes}
        directed_smaller.in_neighbors = {node:set(nx_smaller.predecessors(node)) for node in temp_nodes}
        directed_smaller.out_neighbors = {node:set(nx_smaller.neighbors(node)) for node in temp_nodes}
        # save removable nodes from smaller graph (nodes outside the anchor)
        counter = 0
        for node in temp_nodes:
            if(params.anchor_G.find(node) == params.anchor_G.end()):
                all_removable_nodes[counter] = node
                counter = counter + 1
        # prepare information of bigger directed graph
        directed_bigger.loops = list(nx.nodes_with_selfloops(nx_bigger))
        temp_nodes = list(nx_bigger.nodes())
        directed_bigger.nodes = {node:nx_bigger.nodes[node]["GMNL"] for node in temp_nodes}
        undirected_copy = deepcopy(nx_bigger)
        undirected_copy = undirected_copy.to_undirected()
        directed_bigger.connectivity_neighbors = {node:set(undirected_copy.neighbors(node)) for node in temp_nodes}
        directed_bigger.in_neighbors = {node:set(nx_bigger.predecessors(node)) for node in temp_nodes}
        directed_bigger.out_neighbors = {node:set(nx_bigger.neighbors(node)) for node in temp_nodes}
        # out degrees for ordering
        node_degrees = {node:nx_bigger.out_degree(node) for node in temp_nodes}
        reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_H,
                                                                     temp_nodes,
                                                                     directed_bigger.connectivity_neighbors,
                                                                     node_degrees,
                                                                     all_nodes_bigger)
    else:
        # prepare nodes and neighbors of smaller undirected graph
        undirected_smaller.loops = list(nx.nodes_with_selfloops(nx_smaller))
        temp_nodes = list(nx_smaller.nodes())
        undirected_smaller.nodes = {node:nx_smaller.nodes[node]["GMNL"] for node in temp_nodes}
        undirected_smaller.neighbors = {node:set(nx_smaller.neighbors(node)) for node in temp_nodes}
        # save removable nodes from smaller graph (nodes outside the anchor)
        counter = 0
        for node in temp_nodes:
            if(params.anchor_G.find(node) == params.anchor_G.end()):
                all_removable_nodes[counter] = node
                counter = counter + 1
        # prepare information of bigger undirected graph
        undirected_bigger.loops = list(nx.nodes_with_selfloops(nx_bigger))
        temp_nodes = list(nx_bigger.nodes())
        undirected_bigger.nodes = {node:nx_bigger.nodes[node]["GMNL"] for node in temp_nodes}
        undirected_bigger.neighbors = {node:set(nx_bigger.neighbors(node)) for node in temp_nodes}
        # degrees for ordering
        node_degrees = {node:nx_bigger.degree(node) for node in temp_nodes}
        reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_H,
                                                                     temp_nodes,
                                                                     undirected_bigger.neighbors,
                                                                     node_degrees,
                                                                     all_nodes_bigger)

    # prepare total orders for bigger graph
    # NOTE: inside the intensive routines and in associated search-parameters, the
    # domain graph is called G, while the codomain graph is always called H.
    if(params.directed_graphs):
        counter = 0
        for node in all_nodes_bigger:
            counter = counter + 1
            params.total_order_H[node] = counter
            params.inverse_total_order_H[counter] = node
            params.unmatched_H_backup.insert(node)
            params.unmatched_H_ordered_backup.push_back(node)
    else:
        counter = 0
        for node in all_nodes_bigger:
            counter = counter + 1
            params.total_order_H[node] = counter
            params.inverse_total_order_H[counter] = node
            params.unmatched_H_backup.insert(node)
            params.unmatched_H_ordered_backup.push_back(node)

    # prepare edges
    if(params.edge_labels):
        if(params.directed_graphs):
            # edges for smaller directed graph
            all_edges_smaller = [(node1, node2, info["GMEL"]) for (node1, node2, info) in nx_smaller.edges(data = True)]
            for each_vector in all_edges_smaller:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                directed_smaller.raw_edges.push_back(temp_pair)
                # save labeled edge
                temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                directed_smaller.edges[temp_str] = each_vector[2]
            # edges for bigger directed graph
            all_edges_bigger = [(node1, node2, info["GMEL"]) for (node1, node2, info) in nx_bigger.edges(data = True)]
            for each_vector in all_edges_bigger:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                directed_bigger.raw_edges.push_back(temp_pair)
                # save labeled edge
                temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                directed_bigger.edges[temp_str] = each_vector[2]
        else:
            # edges for smaller undirected graph
            all_edges_smaller = [(node1, node2, info["GMEL"]) for (node1, node2, info) in nx_smaller.edges(data = True)]
            for each_vector in all_edges_smaller:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                undirected_smaller.raw_edges.push_back(temp_pair)
                # save labeled edge
                if(each_vector[0] == each_vector[1]):
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[0])
                    undirected_smaller.edges[temp_str] = each_vector[2]
                else:
                    # save the two label edges to simplify future access
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                    undirected_smaller.edges[temp_str] = each_vector[2]
                    temp_str = to_string(each_vector[1]) + comma + to_string(each_vector[0])
                    undirected_smaller.edges[temp_str] = each_vector[2]
            # edges for bigger undirected graph
            all_edges_bigger = [(node1, node2, info["GMEL"]) for (node1, node2, info) in nx_bigger.edges(data = True)]
            for each_vector in all_edges_bigger:
                # save raw edge
                temp_pair.first = each_vector[0]
                temp_pair.second = each_vector[1]
                undirected_bigger.raw_edges.push_back(temp_pair)
                # save labeled edge
                if(each_vector[0] == each_vector[1]):
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[0])
                    undirected_bigger.edges[temp_str] = each_vector[2]
                else:
                    # save the two label edges to simplify future access
                    temp_str = to_string(each_vector[0]) + comma + to_string(each_vector[1])
                    undirected_bigger.edges[temp_str] = each_vector[2]
                    temp_str = to_string(each_vector[1]) + comma + to_string(each_vector[0])
                    undirected_bigger.edges[temp_str] = each_vector[2]

    # prepare ordered neighbors of bigger graph
    # NOTE: inside the intensive routines and in associated search-parameters, the
    # domain graph is called G, while the codomain graph is always called H.
    if(params.directed_graphs):
        partial_maps_order_neighbors(all_nodes_bigger, directed_bigger.in_neighbors, directed_bigger.in_neighbors_ordered,
                                     params.total_order_H, params.inverse_total_order_H)
        partial_maps_order_neighbors(all_nodes_bigger, directed_bigger.out_neighbors, directed_bigger.out_neighbors_ordered,
                                     params.total_order_H, params.inverse_total_order_H)
    else:
        partial_maps_order_neighbors(all_nodes_bigger, undirected_bigger.neighbors, undirected_bigger.neighbors_ordered,
                                     params.total_order_H, params.inverse_total_order_H)

    # set recursion limit
    required_limit = max(nx_G.order(), nx_H.order())
    current_limit = getrecursionlimit()
    if(current_limit < (scalation_value * required_limit)):
        setrecursionlimit(int(scalation_value * required_limit))

    # obtain maximum reachable extension through iterative trimming
    partial_maps_iterative_trimming(all_removable_nodes,
                                    encoded_extensions, params,
                                    directed_bigger, directed_smaller,
                                    undirected_bigger, undirected_smaller)

    # recover extensions
    if(twisted):
        # untwist encoded extensios
        for each_extension in encoded_extensions:
            # reinitialize extensions holder
            temp_extension.clear()
            # untwist extension
            for each_pair in each_extension:
                temp_pair.first = each_pair.second
                temp_pair.second = each_pair.first
                temp_extension.insert(temp_pair)
            # save untwisted extension
            untwisted_extensions.push_back(temp_extension)
        # decode extensions into python objects
        for each_extension in untwisted_extensions:
            extensions.append(decode_match(list(each_extension), encoded_node_names))
    else:
        # decode extensions python objects
        for each_extension in encoded_extensions:
            extensions.append(decode_match(list(each_extension), encoded_node_names))

    # check if there were any proper extensions, if not return the anchor and the value false
    if(len(extensions) > 0):
        found_proper_extensions = True
    else:
        found_proper_extensions = False
        extensions = [input_anchor]

    # end of function
    return(extensions, found_proper_extensions)





# function: iterative trimming for VF2-MCS-like search -------------------------
cdef void partial_maps_iterative_trimming(cpp_unordered_map[int, int] & all_removable_nodes,
                                          cpp_vector[cpp_set[cpp_pair[int, int]]] & encoded_extensions,
                                          partial_maps_search_params & params,
                                          partial_maps_directed_graph & directed_bigger,
                                          partial_maps_directed_graph & directed_smaller,
                                          partial_maps_undirected_graph & undirected_bigger,
                                          partial_maps_undirected_graph & undirected_smaller) noexcept:

    # local variables
    cdef int i = 0
    cdef int k = 0
    cdef int new_index = 0
    cdef int last_added = 0
    cdef int total_removable = 0
    cdef cpp_vector[int] each_subset
    cdef cpp_vector[int] test_subset
    cdef cpp_vector[cpp_vector[int]] old_subsets
    cdef cpp_vector[cpp_vector[int]] new_subsets

    # prepare data for iteration
    total_removable = <int>(all_removable_nodes.size())

    # test "removal" of empty set (only for graphs with different order)
    if(undirected_smaller.nodes.size() < undirected_bigger.nodes.size()):
        # prepare empty test subset
        test_subset.clear()
        # test smaller graph as induced subgraph of bigger
        if(params.directed_graphs):
            trimm_and_test_subgraph_isomorphism_directed(test_subset,
                                                         all_removable_nodes,
                                                         encoded_extensions,
                                                         params,
                                                         directed_bigger,
                                                         directed_smaller)
        else:
            trimm_and_test_subgraph_isomorphism_undirected(test_subset,
                                                           all_removable_nodes,
                                                           encoded_extensions,
                                                           params,
                                                           undirected_bigger,
                                                           undirected_smaller)
        # finish if extensions were found; next trimmings produce smaller graphs
        if(not encoded_extensions.empty()):
            return

    # test removing non-empty sets only if there are at least 2 removable nodes
    if(total_removable >= 2):

        # test removal of singleton sets
        for i in range(total_removable):
            # créate singleton
            test_subset.clear()
            test_subset.push_back(i)
            # trimm smaller graph and test induced subgraph isomorphism
            if(params.directed_graphs):
                trimm_and_test_subgraph_isomorphism_directed(test_subset,
                                                             all_removable_nodes,
                                                             encoded_extensions,
                                                             params,
                                                             directed_bigger,
                                                             directed_smaller)
            else:
                trimm_and_test_subgraph_isomorphism_undirected(test_subset,
                                                               all_removable_nodes,
                                                               encoded_extensions,
                                                               params,
                                                               undirected_bigger,
                                                               undirected_smaller)
            # finish if one extension was found and no more are required
            if(not encoded_extensions.empty()):
                if(not params.all_extensions):
                    return
            # save singleton
            old_subsets.push_back(test_subset)
        # finish if extensions were found; next trimmings produce smaller graphs
        if(not encoded_extensions.empty()):
            return

        # test removal of sets with more than one element, and up to cardinality
        # N-1 for N removable nodes, since removing the N nodes can only produce
        # the trivial extension, i.e., return the input anchor itself
        for k in range(2, total_removable):
            # clear new subsets holder
            new_subsets.clear()
            # generate new subsets with one more element
            for each_subset in old_subsets:
                last_added = each_subset.back()
                if(last_added < total_removable):
                    for new_index in range(last_added + 1, total_removable):
                        # create subset
                        test_subset = each_subset
                        test_subset.push_back(new_index)
                        # trimm smaller graph and test induced subgraph isomorphism
                        if(params.directed_graphs):
                            trimm_and_test_subgraph_isomorphism_directed(test_subset,
                                                                         all_removable_nodes,
                                                                         encoded_extensions,
                                                                         params,
                                                                         directed_bigger,
                                                                         directed_smaller)
                        else:
                            trimm_and_test_subgraph_isomorphism_undirected(test_subset,
                                                                           all_removable_nodes,
                                                                           encoded_extensions,
                                                                           params,
                                                                           undirected_bigger,
                                                                           undirected_smaller)
                        # finish if one extension was found and no more are required
                        if(not encoded_extensions.empty()):
                            if(not params.all_extensions):
                                return
                        # save subset
                        new_subsets.push_back(test_subset)
            # update subsets holders
            old_subsets.clear()
            old_subsets = new_subsets
            # finish if extensions were found; next trimmings produce smaller graphs
            if(not encoded_extensions.empty()):
                return

    # end of function





# function: prepare smaller graph and test subgraph isomorphism ----------------
cdef void trimm_and_test_subgraph_isomorphism_undirected(cpp_vector[int] & test_subset,
                                                         cpp_unordered_map[int, int] & all_removable_nodes,
                                                         cpp_vector[cpp_set[cpp_pair[int, int]]] & encoded_extensions,
                                                         partial_maps_search_params & params,
                                                         partial_maps_undirected_graph & undirected_bigger,
                                                         partial_maps_undirected_graph & undirected_smaller) noexcept:

    # local variables
    cdef int i = 0
    cdef int deg = 0
    cdef int node = 0
    cdef int offset = 0
    cdef int covered = 0
    cdef int counter = 0
    cdef cpp_bool reachable = True
    cdef cpp_bool consistent_labels = True
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, cpp_unordered_set[int]] each_neighborhood
    cdef cpp_vector[int] deg_trimmed
    cdef cpp_vector[int] nodes_in_trimmed
    cdef cpp_vector[int] ordered_nodes_in_trimmed
    cdef cpp_unordered_set[int] new_neighborhood
    cdef cpp_unordered_set[int] nodes_to_be_removed
    cdef cpp_unordered_map[int, int] node_degrees
    cdef partial_maps_undirected_graph trimmed_graph
    cdef partial_maps_state_undirected initial_state_undirected

    # build trimmed graph only if test subset is not empty
    if(test_subset.empty()):

        # copy graph
        trimmed_graph = undirected_smaller

        # copy degree sequence
        deg_trimmed = undirected_smaller.degrees

        # get nodes and their degree for ordering
        for each_neighborhood in trimmed_graph.neighbors:
            # save node
            node = each_neighborhood.first
            nodes_in_trimmed.push_back(node)
            # save node with its degree
            deg = <int>(each_neighborhood.second.size())
            node_degrees[each_neighborhood.first] = deg

    else:

        # get nodes to be removed
        for i in test_subset:
            nodes_to_be_removed.insert(all_removable_nodes[i])

        # prepare nodes of trimmed graph
        for each_pair in undirected_smaller.nodes:
            # add node only if it is not going to be removed
            if(nodes_to_be_removed.find(each_pair.first) == nodes_to_be_removed.end()):
                trimmed_graph.nodes.insert(each_pair)
                nodes_in_trimmed.push_back(each_pair.first)

        # prepare looped nodes of trimmed graph
        for node in undirected_smaller.loops:
            # add node only if it is not going to be removed
            if(nodes_to_be_removed.find(node) == nodes_to_be_removed.end()):
                trimmed_graph.loops.insert(node)

        # prepare edges nodes of trimmed graph
        for each_pair in undirected_smaller.raw_edges:
            # add edge only if its ends are not nodes to be removed
            if(nodes_to_be_removed.find(each_pair.first) == nodes_to_be_removed.end()):
                if(nodes_to_be_removed.find(each_pair.second) == nodes_to_be_removed.end()):
                    # add raw edge
                    trimmed_graph.raw_edges.push_back(each_pair)
                    # add label edge
                    if(each_pair.first == each_pair.second):
                        temp_str = to_string(each_pair.first) + comma + to_string(each_pair.second)
                        trimmed_graph.edges[temp_str] = undirected_smaller.edges[temp_str]
                    else:
                        temp_str = to_string(each_pair.first) + comma + to_string(each_pair.second)
                        trimmed_graph.edges[temp_str] = undirected_smaller.edges[temp_str]
                        temp_str = to_string(each_pair.second) + comma + to_string(each_pair.first)
                        trimmed_graph.edges[temp_str] = undirected_smaller.edges[temp_str]

        # prepare neighborhoods for trimmed graph
        for each_neighborhood in undirected_smaller.neighbors:
            # add neighborhood only if node is not to be removed
            if(nodes_to_be_removed.find(each_neighborhood.first) == nodes_to_be_removed.end()):
                # initialize new neighborhood
                new_neighborhood.clear()
                # add all nodes that are not to be removed
                for node in each_neighborhood.second:
                    if(nodes_to_be_removed.find(node) == nodes_to_be_removed.end()):
                        new_neighborhood.insert(node)
                # assign neighborhood
                trimmed_graph.neighbors[each_neighborhood.first] = new_neighborhood
                # save degree
                deg = <int>(new_neighborhood.size())
                deg_trimmed.push_back(deg)
                node_degrees[each_neighborhood.first] = deg

    # get order of graphs
    params.expected_order = trimmed_graph.nodes.size()
    params.expected_order_int = <int>(trimmed_graph.nodes.size())
    params.int_order_domain = <int>(trimmed_graph.nodes.size())

    # test for degree cover if necessary
    if(params.test_unbalanced_degree_cover):
        # sort degrees of trimmed graph if necessary
        if(not test_subset.empty()):
            sort(deg_trimmed.begin(), deg_trimmed.end())
        # test if degrees of bigger cover degrees of trimmed
        i = 0
        offset = 0
        covered = 0
        for deg in deg_trimmed:
            # decide if continue searching
            if(offset >= params.int_order_codomain):
                break
            else:
                # search for covering degree
                for i in range(offset, params.int_order_codomain):
                    if(deg <= undirected_bigger.degrees[i]):
                        covered = covered + 1
                        break
                # get next offset
                offset = i + 1
        # test for cover of all degrees
        if(covered < params.expected_order_int):
            return

    # test reachability and concentric order traversal
    reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_G,
                                                                 nodes_in_trimmed,
                                                                 trimmed_graph.neighbors,
                                                                 node_degrees,
                                                                 ordered_nodes_in_trimmed)

    # discard if reachability is requested and it doesnt hold for the trimmed graph
    if(params.reachability):
        if(not reachable):
            return

    # reinitialize parameters of node and edge labes
    params.node_labels = params.node_labels_backup
    params.edge_labels = params.edge_labels_backup

    # test consistency of labels if required
    if(params.node_labels or params.edge_labels):
        consistent_labels = partial_maps_label_consistency(params.caller,
                                                           params.node_labels,
                                                           params.edge_labels,
                                                           trimmed_graph.nodes,
                                                           undirected_bigger.nodes,
                                                           trimmed_graph.edges,
                                                           undirected_bigger.edges,
                                                           trimmed_graph.raw_edges,
                                                           undirected_bigger.raw_edges,
                                                           params.anchor_G,
                                                           params.anchor_H)
        # if labels are not consistent then the reminder graphs cannot be isomorphic
        if(not consistent_labels):
            return

    # reinitialize information in search state for bigger graph
    initial_state_undirected.unmatched_H = params.unmatched_H_backup
    initial_state_undirected.unmatched_H_ordered = params.unmatched_H_ordered_backup

    # prepare total order for nodes in trimmed graph
    params.total_order_G.clear()
    params.inverse_total_order_G.clear()
    counter = 0
    for node in ordered_nodes_in_trimmed:
        counter = counter + 1
        params.total_order_G[node] = counter
        params.inverse_total_order_G[counter] = node
        initial_state_undirected.unmatched_G.insert(node)
        initial_state_undirected.unmatched_G_ordered.push_back(node)

    # prepare ordered neighbors of trimmed graph
    partial_maps_order_neighbors(ordered_nodes_in_trimmed,
                                 trimmed_graph.neighbors,
                                 trimmed_graph.neighbors_ordered,
                                 params.total_order_G,
                                 params.inverse_total_order_G)

    # prepare initial state updated with anchor
    partial_maps_prepare_initial_state_undirected(params.encoded_anchor,
                                                  params.anchor_G,
                                                  params.anchor_H,
                                                  trimmed_graph,
                                                  undirected_bigger,
                                                  initial_state_undirected)

    # evaluate subgraph isomorphism
    partial_maps_undirected(params, initial_state_undirected, trimmed_graph, undirected_bigger, encoded_extensions)

    # end of function





# function: prepare smaller graph and test subgraph isomorphism ----------------
cdef void trimm_and_test_subgraph_isomorphism_directed(cpp_vector[int] & test_subset,
                                                       cpp_unordered_map[int, int] & all_removable_nodes,
                                                       cpp_vector[cpp_set[cpp_pair[int, int]]] & encoded_extensions,
                                                       partial_maps_search_params & params,
                                                       partial_maps_directed_graph & directed_bigger,
                                                       partial_maps_directed_graph & directed_smaller) noexcept:

    # local variables
    cdef int i = 0
    cdef int deg = 0
    cdef int node = 0
    cdef int offset = 0
    cdef int covered = 0
    cdef int counter = 0
    cdef cpp_bool reachable = True
    cdef cpp_bool consistent_labels = True
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, cpp_unordered_set[int]] each_neighborhood
    cdef cpp_vector[int] in_deg_trimmed
    cdef cpp_vector[int] out_deg_trimmed
    cdef cpp_vector[int] nodes_in_trimmed
    cdef cpp_vector[int] ordered_nodes_in_trimmed
    cdef cpp_unordered_set[int] new_neighborhood
    cdef cpp_unordered_set[int] nodes_to_be_removed
    cdef cpp_unordered_map[int, int] node_degrees
    cdef partial_maps_directed_graph trimmed_graph
    cdef partial_maps_state_directed initial_state_directed

    # build trimmed graph only if test subset is not empty
    if(test_subset.empty()):

        # copy graph
        trimmed_graph = directed_smaller

        # copy degree sequences
        in_deg_trimmed = directed_smaller.in_degrees
        out_deg_trimmed = directed_smaller.out_degrees

        # get nodes and their out degree for ordering
        for each_neighborhood in trimmed_graph.out_neighbors:
            # save node
            node = each_neighborhood.first
            nodes_in_trimmed.push_back(node)
            # save node with its degree
            deg = <int>(each_neighborhood.second.size())
            node_degrees[each_neighborhood.first] = deg

    else:

        # get nodes to be removed
        for i in test_subset:
            nodes_to_be_removed.insert(all_removable_nodes[i])

        # prepare nodes of trimmed graph
        for each_pair in directed_smaller.nodes:
            # add node only if it is not going to be removed
            if(nodes_to_be_removed.find(each_pair.first) == nodes_to_be_removed.end()):
                trimmed_graph.nodes.insert(each_pair)
                nodes_in_trimmed.push_back(each_pair.first)

        # prepare looped nodes of trimmed graph
        for node in directed_smaller.loops:
            # add node only if it is not going to be removed
            if(nodes_to_be_removed.find(node) == nodes_to_be_removed.end()):
                trimmed_graph.loops.insert(node)

        # prepare edges nodes of trimmed graph
        for each_pair in directed_smaller.raw_edges:
            # add edge only if its ends are not nodes to be removed
            if(nodes_to_be_removed.find(each_pair.first) == nodes_to_be_removed.end()):
                if(nodes_to_be_removed.find(each_pair.second) == nodes_to_be_removed.end()):
                    # add raw edge
                    trimmed_graph.raw_edges.push_back(each_pair)
                    # add label edge
                    temp_str = to_string(each_pair.first) + comma + to_string(each_pair.second)
                    trimmed_graph.edges[temp_str] = directed_smaller.edges[temp_str]

        # prepare in-neighborhoods for trimmed graph
        for each_neighborhood in directed_smaller.in_neighbors:
            # add in-neighborhood only if node is not to be removed
            if(nodes_to_be_removed.find(each_neighborhood.first) == nodes_to_be_removed.end()):
                # initialize new neighborhood
                new_neighborhood.clear()
                # add all nodes that are not to be removed
                for node in each_neighborhood.second:
                    if(nodes_to_be_removed.find(node) == nodes_to_be_removed.end()):
                        new_neighborhood.insert(node)
                # assign neighborhood
                trimmed_graph.in_neighbors[each_neighborhood.first] = new_neighborhood
                # save degree
                deg = <int>(new_neighborhood.size())
                in_deg_trimmed.push_back(deg)

        # prepare out-neighborhoods for trimmed graph
        for each_neighborhood in directed_smaller.out_neighbors:
            # add out-neighborhood only if node is not to be removed
            if(nodes_to_be_removed.find(each_neighborhood.first) == nodes_to_be_removed.end()):
                # initialize new neighborhood
                new_neighborhood.clear()
                # add all nodes that are not to be removed
                for node in each_neighborhood.second:
                    if(nodes_to_be_removed.find(node) == nodes_to_be_removed.end()):
                        new_neighborhood.insert(node)
                # assign neighborhood
                trimmed_graph.out_neighbors[each_neighborhood.first] = new_neighborhood
                # save degree
                deg = <int>(new_neighborhood.size())
                out_deg_trimmed.push_back(deg)
                node_degrees[each_neighborhood.first] = deg

        # prepare connectivity-neighborhoods for trimmed graph
        for each_neighborhood in directed_smaller.connectivity_neighbors:
            # add connectivity-neighbors only if node is not to be removed
            if(nodes_to_be_removed.find(each_neighborhood.first) == nodes_to_be_removed.end()):
                # initialize new neighborhood
                new_neighborhood.clear()
                # add all nodes that are not to be removed
                for node in each_neighborhood.second:
                    if(nodes_to_be_removed.find(node) == nodes_to_be_removed.end()):
                        new_neighborhood.insert(node)
                # assign neighborhood
                trimmed_graph.connectivity_neighbors[each_neighborhood.first] = new_neighborhood

    # get order of graphs
    params.expected_order = trimmed_graph.nodes.size()
    params.expected_order_int = <int>(trimmed_graph.nodes.size())
    params.int_order_domain = <int>(trimmed_graph.nodes.size())

    # test for in-degree cover if necessary
    if(params.test_unbalanced_in_degree_cover):
        # sort in-degrees of trimmed graph if necessary
        if(not test_subset.empty()):
            sort(in_deg_trimmed.begin(), in_deg_trimmed.end())
        # test if in-degrees of bigger cover in-degrees of trimmed
        i = 0
        offset = 0
        covered = 0
        for deg in in_deg_trimmed:
            # decide if continue searching
            if(offset >= params.int_order_codomain):
                break
            else:
                # search for covering in-degree
                for i in range(offset, params.int_order_codomain):
                    if(deg <= directed_bigger.in_degrees[i]):
                        covered = covered + 1
                        break
                # get next offset
                offset = i + 1
        # test for cover of all degrees
        if(covered < params.expected_order_int):
            return

    # test for out-degree cover if necessary
    if(params.test_unbalanced_out_degree_cover):
        # sort out-degrees of trimmed graph if necessary
        if(not test_subset.empty()):
            sort(out_deg_trimmed.begin(), out_deg_trimmed.end())
        # test if out-degrees of bigger cover out-degrees of trimmed
        i = 0
        offset = 0
        covered = 0
        for deg in out_deg_trimmed:
            # decide if continue searching
            if(offset >= params.int_order_codomain):
                break
            else:
                # search for covering out-degree
                for i in range(offset, params.int_order_codomain):
                    if(deg <= directed_bigger.out_degrees[i]):
                        covered = covered + 1
                        break
                # get next offset
                offset = i + 1
        # test for cover of all degrees
        if(covered < params.expected_order_int):
            return

    # test reachability and concentric order traversal
    reachable = partial_maps_order_nodes_concentric_reachability(params.anchor_G,
                                                                 nodes_in_trimmed,
                                                                 trimmed_graph.connectivity_neighbors,
                                                                 node_degrees,
                                                                 ordered_nodes_in_trimmed)

    # discard if reachability is requested and it doesnt hold for the trimmed graph
    if(params.reachability):
        if(not reachable):
            return

    # reinitialize parameters of node and edge labes
    params.node_labels = params.node_labels_backup
    params.edge_labels = params.edge_labels_backup

    # test consistency of labels if required
    if(params.node_labels or params.edge_labels):
        consistent_labels = partial_maps_label_consistency(params.caller,
                                                           params.node_labels,
                                                           params.edge_labels,
                                                           trimmed_graph.nodes,
                                                           directed_bigger.nodes,
                                                           trimmed_graph.edges,
                                                           directed_bigger.edges,
                                                           trimmed_graph.raw_edges,
                                                           directed_bigger.raw_edges,
                                                           params.anchor_G,
                                                           params.anchor_H)
        # if labels are not consistent then the reminder graphs cannot be isomorphic
        if(not consistent_labels):
            return

    # reinitialize information in search state for bigger graph
    initial_state_directed.unmatched_H = params.unmatched_H_backup
    initial_state_directed.unmatched_H_ordered = params.unmatched_H_ordered_backup

    # prepare total order for nodes in trimmed graph
    params.total_order_G.clear()
    params.inverse_total_order_G.clear()
    counter = 0
    for node in ordered_nodes_in_trimmed:
        counter = counter + 1
        params.total_order_G[node] = counter
        params.inverse_total_order_G[counter] = node
        initial_state_directed.unmatched_G.insert(node)
        initial_state_directed.unmatched_G_ordered.push_back(node)

    # prepare ordered neighbors of trimmed graph
    partial_maps_order_neighbors(ordered_nodes_in_trimmed,
                                 trimmed_graph.in_neighbors,
                                 trimmed_graph.in_neighbors_ordered,
                                 params.total_order_G,
                                 params.inverse_total_order_G)
    partial_maps_order_neighbors(ordered_nodes_in_trimmed,
                                 trimmed_graph.out_neighbors,
                                 trimmed_graph.out_neighbors_ordered,
                                 params.total_order_G,
                                 params.inverse_total_order_G)

    # prepare initial state updated with anchor
    partial_maps_prepare_initial_state_directed(params.encoded_anchor,
                                                params.anchor_G,
                                                params.anchor_H,
                                                trimmed_graph,
                                                directed_bigger,
                                                initial_state_directed)

    # evaluate subgraph isomorphism
    partial_maps_directed(params, initial_state_directed, trimmed_graph, directed_bigger, encoded_extensions)

    # end of function





# functions - input consistency ################################################





# function: test input correctness ---------------------------------------------
cdef cpp_bool partial_maps_input_correctness(nx_G, nx_H, input_anchor, node_labels, edge_labels, all_extensions, reachability, caller):

    # local variables (python)
    test_int = 0
    test_bool = False
    cdef list nodes_G = []
    cdef list nodes_H = []
    cdef list test_list = [0, 0]
    cdef tuple test_tuple = (0, 0)
    test_undir = nx.Graph()
    test_dir = nx.DiGraph()
    each_node = None
    test_entry = None

    # NOTE: caller flag
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs

    # check that first argument is networkx graph or digraph
    if(type(nx_G) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: first argument must be a networkx graph or digraph."))

    # check that second argument is networkx graph or digraph
    if(type(nx_H) not in [type(test_undir), type(test_dir)]):
        raise(ValueError("gmapache: second argument must be a networkx graph or digraph."))

    # check that third argument is a list
    if(type(input_anchor) not in [type(test_list)]):
        raise(ValueError("gmapache: third argument must be a list."))

    # check that fourth argument is boolean
    if(type(node_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument node_labels must be a boolean variable."))

    # check that fifth argument is boolean
    if(type(edge_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument edge_labels must be a boolean variable."))

    # check that sixth argument is boolean
    if(type(all_extensions) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument all_extensions must be a boolean variable."))

    # check that seventh argument is boolean
    if(type(reachability) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument reachability must be a boolean variable."))

    # check that input graphs are of the same type
    if((nx.is_directed(nx_G)) and (not nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))
    if((not nx.is_directed(nx_G)) and (nx.is_directed(nx_H))):
        raise(ValueError("gmapache: input graphs must be both directed or both undirected."))

    # check that input graphs are not null graphs
    if((nx_G.order() == 0) or (nx_H.order() == 0)):
        raise(ValueError("gmapache: input graphs must have at least one node each."))

    # check that input graphs have the same number of nodes
    if(caller == 0):
        if(not nx_G.order() == nx_H.order()):
            return(False)

    # check correctness of input anchor and its entries
    if(len(input_anchor) == 0):
        raise(ValueError("gmapache: third argument must be a non-empty list of 2-tuples."))
    nodes_G = list(nx_G.nodes())
    nodes_H = list(nx_H.nodes())
    for test_entry in input_anchor:
        if(not type(test_entry) in [type(test_tuple)]):
            raise(ValueError("gmapache: all elements in input anchor must be tuples."))
        if(not len(test_entry) == 2):
            raise(ValueError("gmapache: all tuples in input anchor must be of lenght 2."))
        if(test_entry[0] not in nodes_G):
            raise(ValueError("gmapache: the input anchor is matching a node not present in the first graph."))
        if(test_entry[1] not in nodes_H):
            raise(ValueError("gmapache: the input anchor is matching a node not present in the second graph."))

    # check amount of entries in input anchor
    if(not len(list(set([x for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("gmapache: found repeated elements from first graph in anchor; the input anchor must represent an injective map and should not have repeated elements."))
    if(not len(list(set([y for (x, y) in input_anchor]))) == len(input_anchor)):
        raise(ValueError("gmapache: found repeated elements from second graph in anchor; the input anchor must represent an injective map and should not have repeated elements."))

    # check that input anchor does not cover input graphs completely
    if(caller == 0):
        if(len(input_anchor) == nx_G.order()):
            return(False)

    # check that input anchor does not cover any or both input graphs
    if(caller == 1):
        if((len(input_anchor) == nx_G.order()) or (len(input_anchor) == nx_H.order())):
            return(False)

    # end of function
    return(True)





# function: consistency of node or edge labels ---------------------------------
cdef cpp_bool partial_maps_label_consistency(int caller,
                                             cpp_bool & node_labels,
                                             cpp_bool & edge_labels,
                                             cpp_unordered_map[int, int] & nodes_G,
                                             cpp_unordered_map[int, int] & nodes_H,
                                             cpp_unordered_map[cpp_string, int] & edges_G,
                                             cpp_unordered_map[cpp_string, int] & edges_H,
                                             cpp_vector[cpp_pair[int, int]] & raw_edges_G,
                                             cpp_vector[cpp_pair[int, int]] & raw_edges_H,
                                             cpp_unordered_set[int] & anchor_G,
                                             cpp_unordered_set[int] & anchor_H) noexcept:

    # local variables
    cdef int label = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string some_edge
    cdef cpp_pair[int, int] each_pair
    cdef cpp_pair[int, int] node_and_label
    cdef cpp_pair[int, int] label_and_count
    cdef cpp_pair[cpp_string, int] edge_and_label
    cdef cpp_unordered_map[int, int] count_labels_G
    cdef cpp_unordered_map[int, int] count_labels_H

    # NOTE: caller flag
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs

    # consistency of node labels
    if(node_labels):

        # get node label count on G
        for node_and_label in nodes_G:
            # count only if outside the anchor
            if(anchor_G.find(node_and_label.first) == anchor_G.end()):
                # get node label
                label = node_and_label.second
                # count node label
                if(count_labels_G.find(label) != count_labels_G.end()):
                    count_labels_G[label] = count_labels_G[label] + 1
                else:
                    count_labels_G[label] = 1

        # get node label count on H
        for node_and_label in nodes_H:
            # count only if outside the anchor
            if(anchor_H.find(node_and_label.first) == anchor_H.end()):
                # get node label
                label = node_and_label.second
                # count node label
                if(count_labels_H.find(label) != count_labels_H.end()):
                    count_labels_H[label] = count_labels_H[label] + 1
                else:
                    count_labels_H[label] = 1

        # compare node counts in both graphs per label
        if(caller == 0):
            if(count_labels_G.size() != count_labels_H.size()):
                return(False)
            for label_and_count in count_labels_G:
                if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                    return(False)
                else:
                    if(label_and_count.second != count_labels_H[label_and_count.first]):
                        return(False)
        if(caller == 1):
            if(count_labels_G.size() > count_labels_H.size()):
                return(False)
            for label_and_count in count_labels_G:
                if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                    return(False)
                else:
                    if(label_and_count.second > count_labels_H[label_and_count.first]):
                        return(False)

        # if there is only one node label (IN BOTH GRAPHS) then turn off node-label checking, since it is the same label
        if((count_labels_G.size() == 1) and (count_labels_H.size() == 1)):
            node_labels = False

    # consistency of edge labels
    if(edge_labels):

        # clear counts
        count_labels_G.clear()
        count_labels_H.clear()

        # get edge label count on G
        for each_pair in raw_edges_G:
            # unpack nodes
            node1 = each_pair.first
            node2 = each_pair.second
            # count only if outside the anchor
            if((anchor_G.find(node1) == anchor_G.end()) or (anchor_G.find(node2) == anchor_G.end())):
                # get edge
                some_edge = to_string(node1) + comma + to_string(node2)
                # get edge label
                label = edges_G[some_edge]
                # count edge label
                if(count_labels_G.find(label) != count_labels_G.end()):
                    count_labels_G[label] = count_labels_G[label] + 1
                else:
                    count_labels_G[label] = 1

        # get edge label count on H
        for each_pair in raw_edges_H:
            # unpack nodes
            node1 = each_pair.first
            node2 = each_pair.second
            # count only if outside the anchor
            if((anchor_H.find(node1) == anchor_H.end()) or (anchor_H.find(node2) == anchor_H.end())):
                # get edge
                some_edge = to_string(node1) + comma + to_string(node2)
                # get edge label
                label = edges_H[some_edge]
                # count edge label
                if(count_labels_H.find(label) != count_labels_H.end()):
                    count_labels_H[label] = count_labels_H[label] + 1
                else:
                    count_labels_H[label] = 1

        # compare edge counts in both graphs per label
        if(caller == 0):
            if(count_labels_G.size() != count_labels_H.size()):
                return(False)
            for label_and_count in count_labels_G:
                if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                    return(False)
                else:
                    if(label_and_count.second != count_labels_H[label_and_count.first]):
                        return(False)
        if(caller == 1):
            if(count_labels_G.size() > count_labels_H.size()):
                return(False)
            for label_and_count in count_labels_G:
                if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                    return(False)
                else:
                    if(label_and_count.second > count_labels_H[label_and_count.first]):
                        return(False)

        # if there is only one edge label (IN BOTH GRAPHS) then turn off edge-label checking, since it is the same label
        if((count_labels_G.size() == 1) and (count_labels_H.size() == 1)):
            edge_labels = False

    # end of function
    return(True)





# functions - input preparation ################################################





# function: order nodes by decreasing (out) degree -----------------------------
cdef void partial_maps_order_nodes_by_degree(cpp_vector[cpp_pair[int, int]] & nodes_info,
                                             cpp_vector[int] & all_nodes) noexcept:

    # local variables
    cdef int deg = 0
    cdef int node = 0
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] all_degrees
    cdef cpp_vector[int] temp_vector
    cdef cpp_unordered_map[int, cpp_vector[int]] nodes_by_degree

    # partition node information by node degree
    for each_pair in nodes_info:
        # determine if degree was already encountered
        if(nodes_by_degree.find(each_pair.second) != nodes_by_degree.end()):
            # save new pair associated to corresponding degree
            nodes_by_degree[each_pair.second].push_back(each_pair.first)
        else:
            # save new degree and create vector of pairs with node information
            all_degrees.push_back(each_pair.second)
            temp_vector.clear()
            temp_vector.push_back(each_pair.first)
            nodes_by_degree[each_pair.second] = temp_vector

    # sort encountered degrees in decreasing order
    sort(all_degrees.begin(), all_degrees.end())
    reverse(all_degrees.begin(), all_degrees.end())

    # traverse nodes one to save them by decreasing degree
    for deg in all_degrees:
        for node in nodes_by_degree[deg]:
            all_nodes.push_back(node)

    # end of function





# function: order nodes in concentric sets around the anchor -------------------
cdef cpp_bool partial_maps_order_nodes_concentric_reachability(cpp_unordered_set[int] & initial_level,
                                                               cpp_vector[int] & temp_nodes,
                                                               cpp_unordered_map[int, cpp_unordered_set[int]] & connectivity_neighbors,
                                                               cpp_unordered_map[int, int] & node_degrees,
                                                               cpp_vector[int] & all_nodes) noexcept:

    # local variables
    cdef cpp_bool reachable = True
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef cpp_pair[int, int] temp_pair
    cdef cpp_vector[int] next_level
    cdef cpp_vector[int] current_level
    cdef cpp_vector[int] unmatched_level
    cdef cpp_vector[cpp_pair[int, int]] pre_level
    cdef cpp_unordered_map[int, cpp_bool] visited

    # prepare visited-nodes array
    for node in temp_nodes:
        visited[node] = False

    # prepare nodes for current level
    for node in initial_level:
        temp_pair.first = node
        temp_pair.second = node_degrees[node]
        pre_level.push_back(temp_pair)
        visited[node] = True

    # order and save nodes in current level
    partial_maps_order_nodes_by_degree(pre_level, current_level)
    for node in current_level:
        all_nodes.push_back(node)

    # multi-source DFS
    while(not current_level.empty()):
        # reinitialize next_level
        pre_level.clear()
        next_level.clear()
        # iterate getting immediate neighbors not already ordered
        for node1 in current_level:
            for node2 in connectivity_neighbors[node1]:
                # only assign oreder of not visited yet
                if(not visited[node2]):
                    # prepare nodes for next level
                    temp_pair.first = node2
                    temp_pair.second = node_degrees[node2]
                    pre_level.push_back(temp_pair)
                    visited[node2] = True
        # level management
        partial_maps_order_nodes_by_degree(pre_level, next_level)
        for node in next_level:
            all_nodes.push_back(node)
        # update nodes to be ordered
        current_level.clear()
        current_level = next_level

    # test for unreachable nodes
    pre_level.clear()
    unmatched_level.clear()
    for node in temp_nodes:
        if(not visited[node]):
            temp_pair.first = node
            temp_pair.second = node_degrees[node]
            pre_level.push_back(temp_pair)
            visited[node] = True

    # enumerate unreachable nodes
    if(not pre_level.empty()):
        # set reachability to false
        reachable = False
        # order and save unmatched vertices
        partial_maps_order_nodes_by_degree(pre_level, unmatched_level)
        for node in unmatched_level:
            all_nodes.push_back(node)

    # end of function
    return(reachable)





# function: order neighbors according to the total order -----------------------
cdef void partial_maps_order_neighbors(cpp_vector[int] & all_nodes,
                                       cpp_unordered_map[int, cpp_unordered_set[int]] & all_neighbors,
                                       cpp_unordered_map[int, cpp_vector[int]] & all_neighbors_ordered,
                                       cpp_unordered_map[int, int] & total_order,
                                       cpp_unordered_map[int, int] & inverse_total_order) noexcept:

    # local variables
    cdef int node = 0
    cdef int each_order = 0
    cdef int each_neighbor = 0
    cdef cpp_vector[int] temp_vector
    cdef cpp_vector[int] empty_vector
    cdef cpp_pair[int, int] each_pair

    # order neighbors of each node
    for node in all_nodes:
        # get total orders
        temp_vector.clear()
        for each_neighbor in all_neighbors[node]:
            temp_vector.push_back(total_order[each_neighbor])
        # sort total orders increasingly
        sort(temp_vector.begin(), temp_vector.end())
        # get ordered neighbors
        all_neighbors_ordered[node] = empty_vector
        for each_order in temp_vector:
            all_neighbors_ordered[node].push_back(inverse_total_order[each_order])

    # end of function





# function: prepare initial state with input anchor - undirected ---------------
cdef void partial_maps_prepare_initial_state_undirected(cpp_set[cpp_pair[int, int]] & encoded_anchor,
                                                        cpp_unordered_set[int] & anchor_G,
                                                        cpp_unordered_set[int] & anchor_H,
                                                        partial_maps_undirected_graph & G,
                                                        partial_maps_undirected_graph & H,
                                                        partial_maps_state_undirected & initial_state) noexcept:

    # local variables
    cdef cpp_bool inserted
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef cpp_pair[int, int] each_pair

    # extend initial search state with each pair in the anchor
    for each_pair in encoded_anchor:

        # get individual nodes
        node1 = each_pair.first
        node2 = each_pair.second

        # add new pair
        initial_state.match.insert(each_pair)

        # add node to match in G
        initial_state.match_G.insert(node1)

        # add node to match in H
        initial_state.match_H.insert(node2)

        # upgrade forward map
        initial_state.forward_match[node1] = node2

        # upgrade inverse map
        initial_state.inverse_match[node2] = node1

        # remove node from unmatched nodes in G
        initial_state.unmatched_G.erase(node1)
        initial_state.unmatched_G_ordered.remove(node1)

        # remove node from unmatched nodes in H
        initial_state.unmatched_H.erase(node2)
        initial_state.unmatched_H_ordered.remove(node2)

        # add unmatched neighbors to ring in G if not already there
        for node in G.neighbors_ordered[node1]:
            # we only insert things outside the anchor into the match and ring
            if(anchor_G.find(node) == anchor_G.end()):
                inserted = (initial_state.ring_G.insert(node)).second
                if(inserted):
                    initial_state.ring_G_ordered.push_back(node)

        # add unmatched neighbors to ring in H if not already there
        for node in H.neighbors_ordered[node2]:
            # we only insert things outside the anchor into the match and ring
            if(anchor_H.find(node) == anchor_H.end()):
                inserted = (initial_state.ring_H.insert(node)).second
                if(inserted):
                    initial_state.ring_H_ordered.push_back(node)

    # end of function





# function: prepare initial state with input anchor - directed -----------------
cdef void partial_maps_prepare_initial_state_directed(cpp_set[cpp_pair[int, int]] & encoded_anchor,
                                                      cpp_unordered_set[int] & anchor_G,
                                                      cpp_unordered_set[int] & anchor_H,
                                                      partial_maps_directed_graph & G,
                                                      partial_maps_directed_graph & H,
                                                      partial_maps_state_directed & initial_state) noexcept:

    # local variables
    cdef cpp_bool inserted
    cdef int node = 0
    cdef int node1 = 0
    cdef int node2 = 0
    cdef cpp_pair[int, int] each_pair

    # extend initial search state with each pair in the anchor
    for each_pair in encoded_anchor:

        # get individual nodes
        node1 = each_pair.first
        node2 = each_pair.second

        # add new pair
        initial_state.match.insert(each_pair)

        # add node to match in G
        initial_state.match_G.insert(node1)

        # add node to match in H
        initial_state.match_H.insert(node2)

        # upgrade forward map
        initial_state.forward_match[node1] = node2

        # upgrade inverse map
        initial_state.inverse_match[node2] = node1

        # remove node from unmatched nodes in G
        initial_state.unmatched_G.erase(node1)
        initial_state.unmatched_G_ordered.remove(node1)

        # remove node from unmatched nodes in H
        initial_state.unmatched_H.erase(node2)
        initial_state.unmatched_H_ordered.remove(node2)

        # add unmatched in-neighbors to in-ring in G if not already there
        for node in G.in_neighbors_ordered[node1]:
            # we only insert things outside the anchor into the match and ring
            if(anchor_G.find(node) == anchor_G.end()):
                inserted = (initial_state.in_ring_G.insert(node)).second
                if(inserted):
                    initial_state.in_ring_G_ordered.push_back(node)

        # add unmatched out-neighbors to out-ring in G if not already there
        for node in G.out_neighbors_ordered[node1]:
            # we only insert things outside the anchor into the match and ring
            if(anchor_G.find(node) == anchor_G.end()):
                inserted = (initial_state.out_ring_G.insert(node)).second
                if(inserted):
                    initial_state.out_ring_G_ordered.push_back(node)

        # add unmatched in-neighbors to in-ring in H if not already there
        for node in H.in_neighbors_ordered[node2]:
            # we only insert things outside the anchor into the match and ring
            if(anchor_H.find(node) == anchor_H.end()):
                inserted = (initial_state.in_ring_H.insert(node)).second
                if(inserted):
                    initial_state.in_ring_H_ordered.push_back(node)

        # add unmatched out-neighbors to out-ring in H if not already there
        for node in H.out_neighbors_ordered[node2]:
            # we only insert things outside the anchor into the match and ring
            if(anchor_H.find(node) == anchor_H.end()):
                inserted = (initial_state.out_ring_H.insert(node)).second
                if(inserted):
                    initial_state.out_ring_H_ordered.push_back(node)

    # end of function





# functions - search extensions - undirected ###################################





# function: core routine of VF2-like undirected approach -----------------------
cdef void partial_maps_undirected(partial_maps_search_params & params,
                                  partial_maps_state_undirected & current_state,
                                  partial_maps_undirected_graph & G,
                                  partial_maps_undirected_graph & H,
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
    cdef partial_maps_change_in_state_undirected change_in_state

    # save if extensions was reached
    if(current_state.match.size() == params.expected_order):
        all_matches.push_back(current_state.match)
        if(not params.all_extensions):
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

                # evaluate syntactic feasibility
                syntactic_feasibility = syntactic_feasibility_undirected(matchable_node_G,
                                                                         candidate_node_H,
                                                                         params,
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

                # evaluate semantic feasibility
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
                        partial_maps_undirected(params, current_state, G, H, all_matches)

                        # finish if only one extension was requested and it was already found
                        if(not all_matches.empty()):
                            if(not params.all_extensions):
                                return

                        # restore state unordered sets
                        restore_match_undirected(matchable_node_G, candidate_node_H, change_in_state, current_state)

                        # restore state ordered lists
                        current_state.ring_G_ordered = ring_G_ordered_backup
                        current_state.ring_H_ordered = ring_H_ordered_backup
                        current_state.unmatched_G_ordered = unmatched_G_ordered_backup
                        current_state.unmatched_H_ordered = unmatched_H_ordered_backup

    # end of function





# function: get candidate pairs for undirected extension search ----------------
cdef cpp_pair[int, int] candidates_info_undirected(partial_maps_search_params & params,
                                                   partial_maps_state_undirected & current_state) noexcept:

    # output holders
    cdef cpp_pair[int, int] candidates_info

    # local variables
    cdef int node = 0
    cdef int minimum_node = -1
    cdef int minimum_value = 0

    # initialize with impossible minimum (total order ranges from 1 to expected_order)
    minimum_value = 10 + params.expected_order_int

    # build candidates either with nodes in the rings or with unmatched nodes
    if((not current_state.ring_G.empty()) and (not current_state.ring_H.empty())):

        # get node with minimum order from (ordered) ring in G
        for node in current_state.ring_G_ordered:
            if(params.total_order_G[node] < minimum_value):
                minimum_value = params.total_order_G[node]
                minimum_node = node

        # build output pair
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
                                                         partial_maps_state_undirected & current_state) noexcept:

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
cdef partial_maps_change_in_state_undirected extend_match_undirected(int taken_from,
                                                                     int node1,
                                                                     int node2,
                                                                     partial_maps_undirected_graph & G,
                                                                     partial_maps_undirected_graph & H,
                                                                     partial_maps_state_undirected & current_state) noexcept:

    # output holders
    cdef partial_maps_change_in_state_undirected change_in_state

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
                                   partial_maps_change_in_state_undirected & change_in_state,
                                   partial_maps_state_undirected & current_state) noexcept:

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





# function: evaluate syntactic feasability for extension search ----------------
cdef cpp_bool syntactic_feasibility_undirected(int node1,
                                               int node2,
                                               partial_maps_search_params & params,
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

    # NOTE: caller flag
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs

    # consistency of degree
    if(params.caller == 0):
        if(neigh_G[node1].size() != neigh_H[node2].size()):
            return(False)
    if(params.caller == 1):
        if(neigh_G[node1].size() > neigh_H[node2].size()):
            return(False)

    # loop-consistency-test
    if(loops_G.find(node1) != loops_G.end()):
        if(loops_H.find(node2) == loops_H.end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(params.induced_subgraph):
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

    if(params.induced_subgraph):
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
    else:
        for node in ordered_neigh_H[node2]:
            # consistency of match was already checked for non-induced search
            if(current_match_H.find(node) == current_match_H.end()):
                if(ring_H.find(node) != ring_H.end()):
                    # save neighbor since we are just comparing numbers later
                    neighbors_ring_H = neighbors_ring_H + 1

    # look ahead 1: consistency of neighbors in ring (not in match but adjacent to match)
    if(params.caller == 0):
        if(neighbors_ring_G != neighbors_ring_H):
            return(False)
    if(params.caller == 1):
        if(neighbors_ring_G > neighbors_ring_H):
            return(False)

    # look ahead 2: consistency of extern neighbors (neither in match nor adjacent to match)
    # extern neighbors are not preserved in non-induced case, because non-edges are not necessarily preserved
    if(params.induced_subgraph):
        if(params.caller == 0):
            if(neighbors_extern_G != neighbors_extern_H):
                return(False)
        if(params.caller == 1):
            if(neighbors_extern_G > neighbors_extern_H):
                return(False)

    # end of function
    return(True)





# function: evaluate semantic feasability for extension search -----------------
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





# functions - search extensions - directed #####################################





# function: core routine of VF2-like directed approach -------------------------
cdef void partial_maps_directed(partial_maps_search_params & params,
                                partial_maps_state_directed & current_state,
                                partial_maps_directed_graph & G,
                                partial_maps_directed_graph & H,
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
    cdef partial_maps_change_in_state_directed change_in_state

    # save if extension was reached
    if(current_state.match.size() == params.expected_order):
        all_matches.push_back(current_state.match)
        if(not params.all_extensions):
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

                # evaluate in syntactic feasibility
                in_syntactic_feasibility = syntactic_feasibility_directed(matchable_node_G,
                                                                          candidate_node_H,
                                                                          params,
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

                # evaluate out syntactic feasibility
                if(in_syntactic_feasibility):
                    out_syntactic_feasibility = syntactic_feasibility_directed(matchable_node_G,
                                                                               candidate_node_H,
                                                                               params,
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

                    # evaluate semantic feasibility
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
                            partial_maps_directed(params, current_state, G, H, all_matches)

                            # finish if only one extension was requested and it was already found
                            if(not all_matches.empty()):
                                if(not params.all_extensions):
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





# function: get candidate pairs for directed extension search ------------------
cdef cpp_pair[int, int] candidates_info_directed(partial_maps_search_params & params,
                                                 partial_maps_state_directed & current_state) noexcept:

    # output holders
    cdef cpp_pair[int, int] candidates_info

    # local variables
    cdef int node = 0
    cdef int minimum_node = -1
    cdef int minimum_value = 0

    # initialize with impossible minimum (total order ranges from 1 to expected_order)
    minimum_value = 10 + params.expected_order_int

    # build candidates from out ring
    if((not current_state.out_ring_G.empty()) and (not current_state.out_ring_H.empty())):

        # get node with minimum order from out ring in G
        for node in current_state.out_ring_G_ordered:
            if(params.total_order_G[node] < minimum_value):
                minimum_value = params.total_order_G[node]
                minimum_node = node

        # build output pair
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
                                                       partial_maps_state_directed & current_state,
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
cdef partial_maps_change_in_state_directed extend_match_directed(int taken_from,
                                                                 int node1,
                                                                 int node2,
                                                                 partial_maps_directed_graph & G,
                                                                 partial_maps_directed_graph & H,
                                                                 partial_maps_state_directed & current_state) noexcept:

    # output holders
    cdef partial_maps_change_in_state_directed change_in_state

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
                                 partial_maps_change_in_state_directed & change_in_state,
                                 partial_maps_state_directed & current_state) noexcept:

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





# function: evaluate syntactic feasability for extension search ----------------
cdef cpp_bool syntactic_feasibility_directed(int node1,
                                             int node2,
                                             partial_maps_search_params & params,
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

    # NOTE: caller flag
    # caller = 0 -> gm.partial_maps.search_complete_induced_extension
    # caller = 1 -> gm.partial_maps.search_maximum_common_anchored_subgraphs

    # consistency of degrees
    if(params.caller == 0):
        if(neigh_G[node1].size() != neigh_H[node2].size()):
            return(False)
    if(params.caller == 1):
        if(neigh_G[node1].size() > neigh_H[node2].size()):
            return(False)

    # loop-consistency-test
    if(loops_G.find(node1) != loops_G.end()):
        if(loops_H.find(node2) == loops_H.end()):
            # node1 has a loop in G but node2 has no loop in H
            return(False)
    if(params.induced_subgraph):
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

    if(params.induced_subgraph):
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
    else:
        for node in ordered_neigh_H[node2]:
            # consistency of match was already checked for non-induced search
            if(current_match_H.find(node) == current_match_H.end()):
                # neighbor in in-ring
                if(in_ring_H.find(node) != in_ring_H.end()):
                    neighbors_in_ring_H = neighbors_in_ring_H + 1
                # neighbor (possibly also) in out-ring
                if(out_ring_H.find(node) != out_ring_H.end()):
                    neighbors_out_ring_H = neighbors_out_ring_H + 1

    # look ahead 1: consistency of neighbors in rings (not in match but adjacent to match)
    if(params.caller == 0):
        if(neighbors_in_ring_G != neighbors_in_ring_H):
            return(False)
        if(neighbors_out_ring_G != neighbors_out_ring_H):
            return(False)
    if(params.caller == 1):
        if(neighbors_in_ring_G > neighbors_in_ring_H):
            return(False)
        if(neighbors_out_ring_G > neighbors_out_ring_H):
            return(False)

    # look ahead 2: consistency of extern neighbors (neither in match nor adjacent to match)
    # extern neighbors are not preserved in non-induced case, because non-edges are not necessarily preserved
    if(params.induced_subgraph):
        if(params.caller == 0):
            if(neighbors_extern_G != neighbors_extern_H):
                return(False)
        if(params.caller == 1):
            if(neighbors_extern_G > neighbors_extern_H):
                return(False)

    # end of function
    return(True)





# function: evaluate semantic feasability for extension search -----------------
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
