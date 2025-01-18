################################################################################
#                                                                              #
# - GranMapache: GRAphs-and-Networks MAPping Applications                      #
#   with Cython and HEuristics                                                 #
#                                                                              #
# - Module: isomorphisms                                                       #
#                                                                              #
# - Description: analysis of isomorphisms and automorphisms.                   #
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
    # ordered neighbors data (for insertion into ring)
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
    # ordered neighbors data (for insertion into ring)
    cpp_unordered_map[int, cpp_vector[int]] in_neighbors_ordered
    cpp_unordered_map[int, cpp_vector[int]] out_neighbors_ordered





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
    # given total order for G
    cpp_bool got_order_for_G
    # given total order for H
    cpp_bool got_order_for_H
    # return all isomorphisms
    cpp_bool all_isomorphisms
    # return analyzing complement
    cpp_bool complement
    # expected amount of matches
    size_t expected_order
    int expected_order_int
    # total order for VF2-like search
    cpp_unordered_map[int, int] total_order_G
    cpp_unordered_map[int, int] total_order_H
    # inverse total order for VF2-like search
    cpp_unordered_map[int, int] inverse_total_order_G
    cpp_unordered_map[int, int] inverse_total_order_H





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





# algorithms ###################################################################





# functions - search isomorphisms - wrapper ####################################





# function: callable wrapper for searching for isomorphisms --------------------
def search_isomorphisms(nx_G = nx.Graph(),           # can also be a networkx DiGraph
                        nx_H = nx.Graph(),           # can also be a networkx DiGraph
                        node_labels = True,          # consider node labels when evaluating isomorphisms
                        edge_labels = True,          # consider edge labels when evaluating isomorphisms
                        all_isomorphisms = False,    # by default stops when finding one isomorphism (if any)
                        total_order_A = dict(),      # custom total order for the nodes of first graph (nx_G)
                        total_order_B = dict()):     # custom total order for the nodes of second graph (nx_H)

    # description
    """
    > description:
    receives two networkx (di-)graphs G and H both directed or both undirected,
    possibly but not necessarily with either node labels and/or edge labels, with
    the same number of nodes and edges, and a boolean variable indicating if the
    function should finish when finding only one isomorphism from G to H (if any),
    or if it should search for all possible isomorphisms from G to H and return them.
    In addition, the VF2-like search requires a total order for the nodes of one of
    the input graphs. Traditionally this is assigned arbitrarily, but here we follow
    the order used by the networkx interface to enumerate the nodes, which usually
    is the order in which the nodes were added to each graph. Nonetheless a custom
    total order can be provided for either or both graphs with the optional input
    dictionaries total_order_A and total_order_B, which can improve the search if
    these orders encode a bijection almost resembling an isomorphism, obtained,
    for example, from a canonicalization algorithm.

    > input:
    * nx_G - first networkx (di)graph being matched.
    * nx_H - second networkx (di)graph being matched.
    * node_labels - boolean indicating if node labels should be considered for
    the search, which is the default behavior, or if they should be ignored.
    * edge_labels - boolean indicating if edge labels should be considered for
    the search, which is the default behavior, or if they should be ignored.
    * all_isomorphisms - boolean variable indicating if the function should stop
    as soon as one isomorphism is found (if any) -default behavior- or if it
    should search for all possible isomorphisms between the graphs.
    * total_order_A - dictionary mapping all the nodes of nx_G into different
    integers from 1 and up to the order of the graphs, representing the custom
    total order for the nodes of first graph.
    * total_order_B - dictionary mapping all the nodes of nx_H into different
    integers from 1 and up to the order of the graphs, representing the custom
    total order for the nodes of second graph.

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
    * gmapache.isomorphisms.search_isomorphisms_label_consistency
    * gmapache.isomorphisms.search_isomorphisms_order_neighbors
    * gmapache.isomorphisms.search_isomorphisms_undirected
    * gmapache.isomorphisms.search_isomorphisms_directed
    """

    # test input correctness
    search_isomorphisms_input_correctness(nx_G, nx_H, node_labels, edge_labels, all_isomorphisms, total_order_A, total_order_B)

    # output holders
    found_isomorphism = False
    cdef list isomorphisms = []

    # fixed threshold parameters
    cdef float limit_edges = 0.5
    cdef float bigger_graphs = 50
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
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string temp_str
    cdef cpp_pair[int, int] each_pair
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
    cdef cpp_vector[cpp_set[cpp_pair[int, int]]] encoded_isomorphisms
    cdef cpp_set[cpp_pair[int, int]] each_isomorphism
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

    # quick test by comparing the degree sequences of input graphs
    params.directed_graphs = nx.is_directed(nx_G)
    if(params.directed_graphs):
        # get degrees
        in_deg_G = [deg for (node_obj, deg) in list(nx_G.in_degree())]
        in_deg_H = [deg for (node_obj, deg) in list(nx_H.in_degree())]
        out_deg_G = [deg for (node_obj, deg) in list(nx_G.out_degree())]
        out_deg_H = [deg for (node_obj, deg) in list(nx_H.out_degree())]
        # sort degrees
        sort(in_deg_G.begin(), in_deg_G.end())
        sort(in_deg_H.begin(), in_deg_H.end())
        sort(out_deg_G.begin(), out_deg_G.end())
        sort(out_deg_H.begin(), out_deg_H.end())
        # compare degree sequences
        if((not in_deg_G == in_deg_H) or (not out_deg_G == out_deg_H)):
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
    params.got_order_for_G = (len(total_order_A) > 0)
    params.got_order_for_H = (len(total_order_B) > 0)
    if(params.directed_graphs):
        params.complement = (nx_G.size() > ceil(((params.expected_order * (params.expected_order-1)) + params.expected_order) * limit_edges)) and (params.expected_order > bigger_graphs)
    else:
        params.complement = (nx_G.size() > ceil(((params.expected_order * (params.expected_order-1)/2) + params.expected_order) * limit_edges)) and (params.expected_order > bigger_graphs)

    # encode graphs
    encoded_graphs, encoded_node_names, encoded_node_labels, encoded_edge_labels = encode_graphs([nx_G, nx_H])

    # get complement if necessary
    if(params.complement):
        complement_G = nx.complement(encoded_graphs[0])
        complement_H = nx.complement(encoded_graphs[1])

    # prepare nodes, neighbors and total order
    if(params.directed_graphs):
        # nodes of directed G
        all_nodes_G = [(node, info["GMNL"]) for (node, info) in encoded_graphs[0].nodes(data = True)]
        directed_G.loops = list(nx.nodes_with_selfloops(encoded_graphs[0]))
        if(params.got_order_for_G):
            # get basic information
            for each_pair in all_nodes_G:
                # save node
                directed_G.nodes.insert(each_pair)
                # neighbors for G
                directed_G.in_neighbors[each_pair.first] = set(encoded_graphs[0].predecessors(each_pair.first))
                directed_G.out_neighbors[each_pair.first] = set(encoded_graphs[0].neighbors(each_pair.first))
                # neighbors for complement of G
                if(params.complement):
                    directed_G.in_neighbors_complement[each_pair.first] = set(complement_G.predecessors(each_pair.first))
                    directed_G.out_neighbors_complement[each_pair.first] = set(complement_G.neighbors(each_pair.first))
                # save total order for the node from input
                params.total_order_G[each_pair.first] = total_order_A[encoded_node_names[each_pair.first]]
                params.inverse_total_order_G[params.total_order_G[each_pair.first]] = each_pair.first
                initial_state_directed.unmatched_G.insert(each_pair.first)
            # get ordered unmatched nodes for G
            for i in range(1, params.expected_order_int + 1):
                initial_state_directed.unmatched_G_ordered.push_back(params.inverse_total_order_G[i])
        else:
            counter = 0
            for each_pair in all_nodes_G:
                # save node
                directed_G.nodes.insert(each_pair)
                # neighbors for G
                directed_G.in_neighbors[each_pair.first] = set(encoded_graphs[0].predecessors(each_pair.first))
                directed_G.out_neighbors[each_pair.first] = set(encoded_graphs[0].neighbors(each_pair.first))
                # neighbors for complement of G
                if(params.complement):
                    directed_G.in_neighbors_complement[each_pair.first] = set(complement_G.predecessors(each_pair.first))
                    directed_G.out_neighbors_complement[each_pair.first] = set(complement_G.neighbors(each_pair.first))
                # save total order for the node
                counter = counter + 1
                params.total_order_G[each_pair.first] = counter
                params.inverse_total_order_G[counter] = each_pair.first
                initial_state_directed.unmatched_G.insert(each_pair.first)
                initial_state_directed.unmatched_G_ordered.push_back(each_pair.first)
        # nodes of directed H
        all_nodes_H = [(node, info["GMNL"]) for (node, info) in encoded_graphs[1].nodes(data = True)]
        directed_H.loops = list(nx.nodes_with_selfloops(encoded_graphs[1]))
        if(params.got_order_for_H):
            for each_pair in all_nodes_H:
                # save node
                directed_H.nodes.insert(each_pair)
                # neighbors for H
                directed_H.in_neighbors[each_pair.first] = set(encoded_graphs[1].predecessors(each_pair.first))
                directed_H.out_neighbors[each_pair.first] = set(encoded_graphs[1].neighbors(each_pair.first))
                # neighbors for complement of H
                if(params.complement):
                    directed_H.in_neighbors_complement[each_pair.first] = set(complement_H.predecessors(each_pair.first))
                    directed_H.out_neighbors_complement[each_pair.first] = set(complement_H.neighbors(each_pair.first))
                # save total order for the node from input
                params.total_order_H[each_pair.first] = total_order_B[encoded_node_names[each_pair.first]]
                params.inverse_total_order_H[params.total_order_H[each_pair.first]] = each_pair.first
                initial_state_directed.unmatched_H.insert(each_pair.first)
            for i in range(1, params.expected_order_int + 1):
                initial_state_directed.unmatched_H_ordered.push_back(params.inverse_total_order_H[i])
        else:
            counter = 0
            for each_pair in all_nodes_H:
                # save node
                directed_H.nodes.insert(each_pair)
                # neighbors for H
                directed_H.in_neighbors[each_pair.first] = set(encoded_graphs[1].predecessors(each_pair.first))
                directed_H.out_neighbors[each_pair.first] = set(encoded_graphs[1].neighbors(each_pair.first))
                # neighbors for complement of H
                if(params.complement):
                    directed_H.in_neighbors_complement[each_pair.first] = set(complement_H.predecessors(each_pair.first))
                    directed_H.out_neighbors_complement[each_pair.first] = set(complement_H.neighbors(each_pair.first))
                # save total order for the node
                counter = counter + 1
                params.total_order_H[each_pair.first] = counter
                params.inverse_total_order_H[counter] = each_pair.first
                initial_state_directed.unmatched_H.insert(each_pair.first)
                initial_state_directed.unmatched_H_ordered.push_back(each_pair.first)
    else:
        # nodes of undirected G
        all_nodes_G = [(node, info["GMNL"]) for (node, info) in encoded_graphs[0].nodes(data = True)]
        undirected_G.loops = list(nx.nodes_with_selfloops(encoded_graphs[0]))
        if(params.got_order_for_G):
            for each_pair in all_nodes_G:
                # save node
                undirected_G.nodes.insert(each_pair)
                # neighbors for G
                undirected_G.neighbors[each_pair.first] = set(encoded_graphs[0].neighbors(each_pair.first))
                # neighbors for complement of G
                if(params.complement):
                    undirected_G.neighbors_complement[each_pair.first] = set(complement_G.neighbors(each_pair.first))
                # save total order for the node from input
                params.total_order_G[each_pair.first] = total_order_A[encoded_node_names[each_pair.first]]
                params.inverse_total_order_G[params.total_order_G[each_pair.first]] = each_pair.first
                initial_state_undirected.unmatched_G.insert(each_pair.first)
            for i in range(1, params.expected_order_int + 1):
                initial_state_undirected.unmatched_G_ordered.push_back(params.inverse_total_order_G[i])
        else:
            counter = 0
            for each_pair in all_nodes_G:
                # save node
                undirected_G.nodes.insert(each_pair)
                # neighbors for G
                undirected_G.neighbors[each_pair.first] = set(encoded_graphs[0].neighbors(each_pair.first))
                # neighbors for complement of G
                if(params.complement):
                    undirected_G.neighbors_complement[each_pair.first] = set(complement_G.neighbors(each_pair.first))
                # save total order for the node
                counter = counter + 1
                params.total_order_G[each_pair.first] = counter
                params.inverse_total_order_G[counter] = each_pair.first
                initial_state_undirected.unmatched_G.insert(each_pair.first)
                initial_state_undirected.unmatched_G_ordered.push_back(each_pair.first)
        # nodes of undirected H
        all_nodes_H = [(node, info["GMNL"]) for (node, info) in encoded_graphs[1].nodes(data = True)]
        undirected_H.loops = list(nx.nodes_with_selfloops(encoded_graphs[1]))
        if(params.got_order_for_H):
            for each_pair in all_nodes_H:
                # save node
                undirected_H.nodes.insert(each_pair)
                # neighbors for H
                undirected_H.neighbors[each_pair.first] = set(encoded_graphs[1].neighbors(each_pair.first))
                # neighbors for complement of H
                if(params.complement):
                    undirected_H.neighbors_complement[each_pair.first] = set(complement_H.neighbors(each_pair.first))
                # save total order for the node from input
                params.total_order_H[each_pair.first] = total_order_B[encoded_node_names[each_pair.first]]
                params.inverse_total_order_H[params.total_order_H[each_pair.first]] = each_pair.first
                initial_state_undirected.unmatched_H.insert(each_pair.first)
            for i in range(1, params.expected_order_int + 1):
                initial_state_undirected.unmatched_H_ordered.push_back(params.inverse_total_order_H[i])
        else:
            counter = 0
            for each_pair in all_nodes_H:
                # save node
                undirected_H.nodes.insert(each_pair)
                # neighbors for H
                undirected_H.neighbors[each_pair.first] = set(encoded_graphs[1].neighbors(each_pair.first))
                # neighbors for complement of H
                if(params.complement):
                    undirected_H.neighbors_complement[each_pair.first] = set(complement_H.neighbors(each_pair.first))
                # save total order for the node
                counter = counter + 1
                params.total_order_H[each_pair.first] = counter
                params.inverse_total_order_H[counter] = each_pair.first
                initial_state_undirected.unmatched_H.insert(each_pair.first)
                initial_state_undirected.unmatched_H_ordered.push_back(each_pair.first)

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

    # test consistency of labels if required
    if(params.node_labels or params.edge_labels):
        # test consistency of node and/or edge labels
        if(params.directed_graphs):
            search_isomorphisms_label_consistency(params.node_labels, params.edge_labels,
                                                  directed_G.nodes, directed_H.nodes,
                                                  directed_G.edges, directed_H.edges)
        else:
            search_isomorphisms_label_consistency(params.node_labels, params.edge_labels,
                                                  undirected_G.nodes, undirected_H.nodes,
                                                  undirected_G.edges, undirected_H.edges)

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

    # evaluate isomorphisms
    if(params.directed_graphs):
        raise(ValueError("gmapache: WORK IN PROGRESS - NOT YET IMPLEMENTED."))
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
cdef void search_isomorphisms_input_correctness(nx_G, nx_H, node_labels, edge_labels, all_isomorphisms, total_order_A, total_order_B):

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

    # check that third argument is networkx graph or digraph
    if(type(node_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument node_labels must be a boolean variable."))

    # check that fourth argument is networkx graph or digraph
    if(type(edge_labels) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument edge_labels must be a boolean variable."))

    # check that fifth argument is networkx graph or digraph
    if(type(all_isomorphisms) not in [type(test_bool)]):
        raise(ValueError("gmapache: argument all_isomorphisms must be a boolean variable."))

    # check that sixth argument is a dictionary
    if(type(total_order_A) not in [type(test_dict)]):
        raise(ValueError("gmapache: argument total_order_A must be a dictionary."))

    # check that seventh argument is a dictionary
    if(type(total_order_B) not in [type(test_dict)]):
        raise(ValueError("gmapache: argument total_order_B must be a dictionary."))

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
        raise(ValueError("gmapache: input graphs must have the same number of nodes."))

    # check that input graphs have the same number of edges
    if(not nx_G.size() == nx_H.size()):
        raise(ValueError("gmapache: input graphs must have the same number of edges."))

    # check consistency of sixth argument
    if(len(total_order_A) > 0):
        if(len(total_order_A) < nx_G.order()):
            raise(ValueError("gmapache: argument total_order_A is missing the map of some nodes from the first graph."))
        if(len(total_order_A) > nx_G.order()):
            raise(ValueError("gmapache: argument total_order_A is mapping more elements than nodes in the first graph."))
        nodes_G = list(nx_G.nodes())
        all_keys = list(total_order_A.keys())
        minimum_order = 10 + len(nodes_G)
        maximum_order = 0
        all_orders.clear()
        for each_node in nodes_G:
            if(each_node not in all_keys):
                raise(ValueError("gmapache: argument total_order_A is missing the map of some nodes from the first graph."))
            if(type(total_order_A[each_node]) not in [type(test_int)]):
                raise(ValueError("gmapache: argument total_order_A must map all the nodes of the first graph to integers."))
            else:
                inserted = all_orders.insert(total_order_A[each_node]).second
                if(not inserted):
                    raise(ValueError("gmapache: argument total_order_A has repeated maps for some of the nodes from the first graph."))
                if(total_order_A[each_node] < minimum_order):
                    minimum_order = total_order_A[each_node]
                if(total_order_A[each_node] > maximum_order):
                    maximum_order = total_order_A[each_node]
        if(minimum_order != 1):
            raise(ValueError("gmapache: argument total_order_A is not a valid total order for the nodes of the first graph."))
        if(maximum_order != len(nodes_G)):
            raise(ValueError("gmapache: argument total_order_A is not a valid total order for the nodes of the first graph."))

    # check consistency of seventh argument
    if(len(total_order_B) > 0):
        if(len(total_order_B) < nx_H.order()):
            raise(ValueError("gmapache: argument total_order_B is missing the map of some nodes from the second graph."))
        if(len(total_order_B) > nx_H.order()):
            raise(ValueError("gmapache: argument total_order_B is mapping more elements than nodes in the second graph."))
        nodes_H = list(nx_H.nodes())
        all_keys = list(total_order_B.keys())
        minimum_order = 10 + len(nodes_H)
        maximum_order = 0
        all_orders.clear()
        for each_node in nodes_H:
            if(each_node not in all_keys):
                raise(ValueError("gmapache: argument total_order_B is missing the map of some nodes from the second graph."))
            if(type(total_order_B[each_node]) not in [type(test_int)]):
                raise(ValueError("gmapache: argument total_order_B must map all the nodes of the second graph to integers."))
            else:
                inserted = all_orders.insert(total_order_B[each_node]).second
                if(not inserted):
                    raise(ValueError("gmapache: argument total_order_B has repeated maps for some of the nodes from the second graph."))
                if(total_order_B[each_node] < minimum_order):
                    minimum_order = total_order_B[each_node]
                if(total_order_B[each_node] > maximum_order):
                    maximum_order = total_order_B[each_node]
        if(minimum_order != 1):
            raise(ValueError("gmapache: argument total_order_B is not a valid total order for the nodes of the second graph."))
        if(maximum_order != len(nodes_H)):
            raise(ValueError("gmapache: argument total_order_B is not a valid total order for the nodes of the second graph."))

    # end of function





# function: consistency of node or edge labels ---------------------------------
cdef void search_isomorphisms_label_consistency(cpp_bool & node_labels,
                                                cpp_bool & edge_labels,
                                                cpp_unordered_map[int, int] & nodes_G,
                                                cpp_unordered_map[int, int] & nodes_H,
                                                cpp_unordered_map[cpp_string, int] & edges_G,
                                                cpp_unordered_map[cpp_string, int] & edges_H):

    # local variables
    cdef int label = 0
    cdef cpp_pair[int, int] node_and_label
    cdef cpp_pair[int, int] label_and_count
    cdef cpp_pair[cpp_string, int] edge_and_label
    cdef cpp_unordered_map[int, int] count_labels_G
    cdef cpp_unordered_map[int, int] count_labels_H

    # consistency of node labels
    if(node_labels):

        # get node label count on G
        for node_and_label in nodes_G:
            # get node label
            label = node_and_label.second
            # count node label
            if(count_labels_G.find(label) != count_labels_G.end()):
                count_labels_G[label] = count_labels_G[label] + 1
            else:
                count_labels_G[label] = 1

        # get node label count on H
        for node_and_label in nodes_H:
            # get node label
            label = node_and_label.second
            # count node label
            if(count_labels_H.find(label) != count_labels_H.end()):
                count_labels_H[label] = count_labels_H[label] + 1
            else:
                count_labels_H[label] = 1

        # compare node counts in both graphs per label
        if(count_labels_G.size() != count_labels_H.size()):
            raise(ValueError("gmapache: requested preservation of node labels but input graphs have different sets of node labels."))
        for label_and_count in count_labels_G:
            if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                raise(ValueError("gmapache: requested preservation of node labels but input graphs have different sets of node labels."))
            else:
                if(label_and_count.second != count_labels_H[label_and_count.first]):
                    raise(ValueError("gmapache: requested preservation of node labels but input graphs have different amount of nodes for some labels."))

        # if there is only one node label then turn off node-label checking
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
            raise(ValueError("gmapache: requested preservation of edge labels but input graphs have different sets of edge labels."))
        for label_and_count in count_labels_G:
            if(count_labels_H.find(label_and_count.first) == count_labels_H.end()):
                raise(ValueError("gmapache: requested preservation of edge labels but input graphs have different sets of edge labels."))
            else:
                if(label_and_count.second != count_labels_H[label_and_count.first]):
                    raise(ValueError("gmapache: requested preservation of edge labels but input graphs have different amount of edges for some labels."))

        # if there is only one edge label then turn off edge-label checking
        if(count_labels_G.size() == 1):
            edge_labels = False

    # end of function





# functions - input preparation ################################################





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

    # order nodes
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





# functions - search isomorphisms - undirected #################################





# function: core routine of VF2-like undirected approach - recursive -----------
cdef void search_isomorphisms_undirected(isomorphisms_search_params & params,
                                         isomorphisms_state_undirected & current_state,
                                         isomorphisms_undirected_graph & G,
                                         isomorphisms_undirected_graph & H,
                                         cpp_vector[cpp_set[cpp_pair[int, int]]] & all_matches) noexcept:

    # local variables
    cdef int node = 0
    cdef int viable_node_G = 0
    cdef int target_node_H = 0
    cdef cpp_bool semantic_feasibility_res = True
    cdef cpp_bool syntactic_feasibility_res = True
    cdef cpp_pair[int, cpp_bool] candidates_info
    cdef cpp_list[int] ordered_candidates
    cdef cpp_list[int] ring_G_ordered_backup
    cdef cpp_list[int] ring_H_ordered_backup
    cdef cpp_list[int] unmatched_G_ordered_backup
    cdef cpp_list[int] unmatched_H_ordered_backup
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

        # get minimum node in H
        target_node_H = candidates_info.first

        # get list of candidates in H
        if(candidates_info.second):
            # from ring
            ordered_candidates = current_state.ring_G_ordered
        else:
            # from unmatched
            ordered_candidates = current_state.unmatched_G_ordered

        # evaluate candidate pairs
        for viable_node_G in ordered_candidates:

            # evaluate syntactic feasibility (possibly over complement graph)
            if(params.complement):
                syntactic_feasibility_res = syntactic_feasibility(viable_node_G,
                                                                  target_node_H,
                                                                  current_state.ring_G,
                                                                  current_state.ring_H,
                                                                  current_state.match_G,
                                                                  current_state.match_H,
                                                                  G.loops,
                                                                  H.loops,
                                                                  current_state.forward_match,
                                                                  current_state.inverse_match,
                                                                  G.neighbors_complement,
                                                                  H.neighbors_complement)
            else:
                syntactic_feasibility_res = syntactic_feasibility(viable_node_G,
                                                                  target_node_H,
                                                                  current_state.ring_G,
                                                                  current_state.ring_H,
                                                                  current_state.match_G,
                                                                  current_state.match_H,
                                                                  G.loops,
                                                                  H.loops,
                                                                  current_state.forward_match,
                                                                  current_state.inverse_match,
                                                                  G.neighbors,
                                                                  H.neighbors)

            # evaluate semantic feasibility (always over original graphs)
            if(syntactic_feasibility_res):
                if(params.node_labels or params.edge_labels):
                    semantic_feasibility_res = semantic_feasibility(params.node_labels,
                                                                    params.edge_labels,
                                                                    viable_node_G,
                                                                    target_node_H,
                                                                    current_state.ring_G,
                                                                    current_state.ring_H,
                                                                    current_state.match_G,
                                                                    current_state.match_H,
                                                                    G.loops,
                                                                    current_state.forward_match,
                                                                    G.nodes,
                                                                    H.nodes,
                                                                    G.neighbors,
                                                                    H.neighbors,
                                                                    G.edges,
                                                                    H.edges)

                # push to stack if valid
                if(semantic_feasibility_res):
                    # extend match with the new pair
                    change_in_state = extend_match_undirected(viable_node_G, target_node_H, G, H, current_state)
                    # recursive call
                    search_isomorphisms_undirected(params, current_state, G, H, all_matches)
                    # finish if only one isosmorphism was requested and it was already found
                    if(not all_matches.empty()):
                        if(not params.all_isomorphisms):
                            return
                    # restore state unordered sets
                    restore_match_undirected(viable_node_G, target_node_H, change_in_state, current_state)
                    # restore state ordered lists
                    current_state.ring_G_ordered = ring_G_ordered_backup
                    current_state.ring_H_ordered = ring_H_ordered_backup
                    current_state.unmatched_G_ordered = unmatched_G_ordered_backup
                    current_state.unmatched_H_ordered = unmatched_H_ordered_backup

    # end of function





# function: get candidate pairs for undirected isomorphism search --------------
cdef cpp_pair[int, cpp_bool] candidates_info_undirected(isomorphisms_search_params & params,
                                                        isomorphisms_state_undirected & current_state) noexcept:

    # output holders
    cdef cpp_pair[int, cpp_bool] candidates_info

    # local variables
    cdef int node = 0
    cdef int minimum_node = -1
    cdef int minimum_value = 0

    # initialize with impossible minimum (total order ranges from 1 to expected_order)
    minimum_value = 10 + params.expected_order_int

    # build candidates either with nodes in the rings or with unmatched nodes
    if((not current_state.ring_G.empty()) and (not current_state.ring_H.empty())):

        # get node with minimum order from ring in H
        for node in current_state.ring_H_ordered:
            if(params.total_order_H[node] < minimum_value):
                minimum_value = params.total_order_H[node]
                minimum_node = node

        # build output pair
        candidates_info.first = minimum_node
        candidates_info.second =  True

    else:

        # get node with minimum order from unmatched nodes in H
        minimum_node = current_state.unmatched_H_ordered.front()

        # build output pair
        candidates_info.first = minimum_node
        candidates_info.second = False

    # end of function
    return(candidates_info)





# function: upgrade search state by extending match ----------------------------
cdef isomorphisms_change_in_state_undirected extend_match_undirected(int node1,
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
    if((not current_state.ring_G.empty()) and (not current_state.ring_H.empty())):
        # remove node from ring in G
        current_state.ring_G_ordered.remove(node1)
        change_in_state.removed_node_ring_G = current_state.ring_G.erase(node1)
        # remove node from ring in H
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

    # add neighbors to ring in G if not already there
    while(not change_in_state.added_neighbors_ring_G.empty()):
        current_state.ring_G.erase(change_in_state.added_neighbors_ring_G.top())
        change_in_state.added_neighbors_ring_G.pop()

    # add neighbors to ring in H if not already there
    while(not change_in_state.added_neighbors_ring_H.empty()):
        current_state.ring_H.erase(change_in_state.added_neighbors_ring_H.top())
        change_in_state.added_neighbors_ring_H.pop()

    # end of function





# functions - search isomorphisms - directed ###################################





# functions - feasability of matches - undirected and directed #################





# function: evaluate syntactic feasability for isomorphism search --------------
cdef cpp_bool syntactic_feasibility(int node1,
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
                                    cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H) noexcept:

    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef int neighbors_ring_G = 0
    cdef int neighbors_ring_H = 0
    cdef int neighbors_extern_G = 0
    cdef int neighbors_extern_H = 0

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
    for node in neigh_G[node1]:
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

    for node in neigh_H[node2]:
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
cdef cpp_bool semantic_feasibility(cpp_bool node_labels,
                                   cpp_bool edge_labels,
                                   int node1,
                                   int node2,
                                   cpp_unordered_set[int] & ring_G,
                                   cpp_unordered_set[int] & ring_H,
                                   cpp_unordered_set[int] & current_match_G,
                                   cpp_unordered_set[int] & current_match_H,
                                   cpp_unordered_set[int] & loops_G,
                                   cpp_unordered_map[int, int] & forward_match,
                                   cpp_unordered_map[int, int] & nodes_G,
                                   cpp_unordered_map[int, int] & nodes_H,
                                   cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_G,
                                   cpp_unordered_map[int, cpp_unordered_set[int]] & neigh_H,
                                   cpp_unordered_map[cpp_string, int] & edges_G,
                                   cpp_unordered_map[cpp_string, int] & edges_H) noexcept:

    # local variables
    cdef int node = 0
    cdef int mapped = 0
    cdef int node_label = 0
    cdef int edge_label = 0
    cdef cpp_string comma
    comma.push_back(44)
    cdef cpp_string labeled_edge_G
    cdef cpp_string labeled_edge_H
    cdef cpp_pair[int, int] each_pair
    cdef cpp_vector[int] neighbors_ring_G
    cdef cpp_vector[int] neighbors_ring_H
    cdef cpp_vector[int] neighbors_match_G
    cdef cpp_vector[int] neighbors_match_H
    cdef cpp_vector[int] neighbors_extern_G
    cdef cpp_vector[int] neighbors_extern_H
    cdef cpp_unordered_map[int, int] count_node_ring_G
    cdef cpp_unordered_map[int, int] count_node_ring_H
    cdef cpp_unordered_map[int, int] count_edge_ring_G
    cdef cpp_unordered_map[int, int] count_edge_ring_H
    cdef cpp_unordered_map[int, int] count_node_extern_G
    cdef cpp_unordered_map[int, int] count_node_extern_H
    cdef cpp_unordered_map[int, int] count_edge_extern_G
    cdef cpp_unordered_map[int, int] count_edge_extern_H

    # compare node-labels
    if(node_labels):
        if(nodes_G[node1] != nodes_H[node2]):
            return(False)

    # compare loop-labels
    if(edge_labels):
        if(loops_G.find(node1) != loops_G.end()):
            labeled_edge_G = to_string(node1) + comma + to_string(node1)
            labeled_edge_H = to_string(node2) + comma + to_string(node2)
            # compare labeled edges
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # obtain tripartition of neighbors in G
    for node in neigh_G[node1]:
        if(current_match_G.find(node) != current_match_G.end()):
            neighbors_match_G.push_back(node)
        else:
            if(ring_G.find(node) != ring_G.end()):
                # save neighbor since we are just comparing numbers later
                neighbors_ring_G.push_back(node)
            else:
                # save neighbor since we are just comparing numbers later
                neighbors_extern_G.push_back(node)

    # label look ahead 0: compare non-loop edge-labels in match
    if(edge_labels):
        for node in neighbors_match_G:
            # edge in G with only one end in match
            labeled_edge_G = to_string(node1) + comma + to_string(node)
            # edge in H with only one end in match
            labeled_edge_H = to_string(node2) + comma + to_string(forward_match[node])
            # compare labeled edges
            if(edges_G[labeled_edge_G] != edges_H[labeled_edge_H]):
                return(False)

    # obtain tripartition of neighbors in H
    for node in neigh_H[node2]:
        if(current_match_H.find(node) != current_match_H.end()):
            neighbors_match_H.push_back(node)
        else:
            if(ring_H.find(node) != ring_H.end()):
                neighbors_ring_H.push_back(node)
            else:
                neighbors_extern_H.push_back(node)

    # label look ahead 1: compare labels of neighbors in ring (not in match but adjacent to match)
    if(not neighbors_ring_G.empty()):
        # count in G
        for node in neighbors_ring_G:
            if(node_labels):
                # get node label
                node_label = nodes_G[node]
                # count node label
                if(count_node_ring_G.find(node_label) != count_node_ring_G.end()):
                    count_node_ring_G[node_label] = count_node_ring_G[node_label] + 1
                else:
                    count_node_ring_G[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                edge_label = edges_G[labeled_edge_G]
                # count edge label
                if(count_edge_ring_G.find(edge_label) != count_edge_ring_G.end()):
                    count_edge_ring_G[edge_label] = count_edge_ring_G[edge_label] + 1
                else:
                    count_edge_ring_G[edge_label] = 1

        # count in H
        for node in neighbors_ring_H:
            if(node_labels):
                # get node label
                node_label = nodes_H[node]
                # count node label
                if(count_node_ring_H.find(node_label) != count_node_ring_H.end()):
                    count_node_ring_H[node_label] = count_node_ring_H[node_label] + 1
                else:
                    count_node_ring_H[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_H = to_string(node2) + comma + to_string(node)
                edge_label = edges_H[labeled_edge_H]
                # count edge label
                if(count_edge_ring_H.find(edge_label) != count_edge_ring_H.end()):
                    count_edge_ring_H[edge_label] = count_edge_ring_H[edge_label] + 1
                else:
                    count_edge_ring_H[edge_label] = 1

        # compare number of types of adjacent nodes
        if(node_labels):
            if(count_node_ring_G.size() != count_node_ring_H.size()):
                return(False)

        # compare number of types of incident edges
        if(edge_labels):
            if(count_edge_ring_G.size() != count_edge_ring_H.size()):
                return(False)

        # compare types of adjacent nodes
        if(node_labels):
            for each_pair in count_node_ring_G:
                if(count_node_ring_H.find(each_pair.first) == count_node_ring_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_node_ring_H[each_pair.first]):
                        return(False)

        # compare types of incident edges
        if(edge_labels):
            for each_pair in count_edge_ring_G:
                if(count_edge_ring_H.find(each_pair.first) == count_edge_ring_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_edge_ring_H[each_pair.first]):
                        return(False)

    # label look ahead 2: compare labels of extern neighbors (neither in match nor adjacent to match)
    if(not neighbors_extern_G.empty()):
        for node in neighbors_extern_G:
            if(node_labels):
                # get node label
                node_label = nodes_G[node]
                # count node label
                if(count_node_extern_G.find(node_label) != count_node_extern_G.end()):
                    count_node_extern_G[node_label] = count_node_extern_G[node_label] + 1
                else:
                    count_node_extern_G[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_G = to_string(node1) + comma + to_string(node)
                edge_label = edges_G[labeled_edge_G]
                # count edge label
                if(count_edge_extern_G.find(edge_label) != count_edge_extern_G.end()):
                    count_edge_extern_G[edge_label] = count_edge_extern_G[edge_label] + 1
                else:
                    count_edge_extern_G[edge_label] = 1

        for node in neighbors_extern_H:
            if(node_labels):
                # get node label
                node_label = nodes_H[node]
                # count node label
                if(count_node_extern_H.find(node_label) != count_node_extern_H.end()):
                    count_node_extern_H[node_label] = count_node_extern_H[node_label] + 1
                else:
                    count_node_extern_H[node_label] = 1
            if(edge_labels):
                # get edge label (possibly loop label)
                labeled_edge_H = to_string(node2) + comma + to_string(node)
                edge_label = edges_H[labeled_edge_H]
                # count edge label
                if(count_edge_extern_H.find(edge_label) != count_edge_extern_H.end()):
                    count_edge_extern_H[edge_label] = count_edge_extern_H[edge_label] + 1
                else:
                    count_edge_extern_H[edge_label] = 1

        # compare number of types of adjacent nodes
        if(node_labels):
            if(count_node_extern_G.size() != count_node_extern_H.size()):
                return(False)

        # compare number of types of incident edges
        if(edge_labels):
            if(count_edge_extern_G.size() != count_edge_extern_H.size()):
                return(False)

        # compare types of adjacent nodes
        if(node_labels):
            for each_pair in count_node_extern_G:
                if(count_node_extern_H.find(each_pair.first) == count_node_extern_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_node_extern_H[each_pair.first]):
                        return(False)

        # compare types of incident edges
        if(edge_labels):
            for each_pair in count_edge_extern_G:
                if(count_edge_extern_H.find(each_pair.first) == count_edge_extern_H.end()):
                    return(False)
                else:
                    if(each_pair.second != count_edge_extern_H[each_pair.first]):
                        return(False)

    # end of function
    return(True)





################################################################################
################################################################################
